# Umsetzungsplan: `apscale/j_generate_read_table.py`

Branch: `refactor/read-table-duckdb`
Arbeitsweise: ein Schritt nach dem anderen. Jeder Schritt = eigener Commit und erst abhaken, wenn
**Test grün** und (wo angegeben) **Benchmark notiert** ist.

Legende je Schritt: `[ ] Impl` · `[ ] Test` · `[ ] Bench` · `[ ] Commit`

---

## Leitplanken

- **Gruppierung muss stabil bleiben.** Gleiche Eingabe → gleiche `sequence_order` → gleiches `group_mapping`, bei jedem Lauf, jeder Threadzahl, jeder Dateireihenfolge.
- **Memory-safe.** Keine Schritte, die Daten proportional zu (Sequenzen × Samples) in Python/pandas halten. Erlaubt: O(Sequenzen) kleine Arrays (z. B. ein int32-Array, 1 Mio. Sequenzen ≈ 4 MB). Alles Große bleibt in DuckDB (spillt auf Platte).
- Öffentliche Funktionssignaturen und Dateinamen der Ausgaben bleiben gleich (`__main__.py` und die Streamlit-Seite rufen `main()` auf).
- Vergleiche in Tests immer **inhaltlich** (Mengen/Sortierung explizit), nicht über Zufallsreihenfolgen.

---

## Phase 0 – Fundament (vor jeder Änderung am Code)

### 0.1 Test- und Benchmark-Gerüst
- [ ] `tests/` anlegen. **Primär echte Daten** (Ordner vom Nutzer, wird nur kopiert, nie verändert); synthetische Daten nur noch für gezielte Randfälle (Ties, Sonderzeichen, leere Samples, >16 384 Samples), die in echten Daten fehlen
- [x] `tests/harness.py`: baut Projekt aus Datensatz-Ordner, führt `main()` aus, misst Laufzeit **und Peak-RSS** (Prozessbaum), schreibt kanonischen Snapshot (nach Hash/Sample, ohne Indizes); `tests/compare_runs.py` vergleicht zwei Snapshots
- [x] **Baseline von `main`** (Commit `d6bb6fa`) in `D:pscale_test_dataaseline_main\{small,small_rerun,large,large_rerun}` (siehe Benchmark-Protokoll unten)
- [ ] Demos aus der Review (Bug 2.1 Gruppen-Parquet) als Test in `tests/` übernehmen (auf echten Daten: `snapshot` vs. Parquet-Zeilenzahl)
- [x] Commit

### 0.2 Reale Daten
- [x] `D:pscale_test_data\small_dataset` (400 Samples, 4 leer, 456 Seqs, 120 Gruppen) und `large_dataset` (639 Samples, 27 964 Seqs, 6 485 Gruppen). Header `hash;size=N` → Pipeline-Stand „nach 10_nc_removal“, kein eigenes Hashing

---

## Phase 1 – Determinismus (Voraussetzung für stabile Gruppierung)

### 1.1 Tie-Breaker + deterministische Indizes (Review 2.4)
Problem (reproduziert): bei gleicher Read-Summe ist `sequence_order` pro Lauf zufällig (im Extremtest 200 000 / 200 000 Sequenzen anders). Betrifft vor allem Singletons, die in Amplikondaten sehr häufig sind.
Lösung: `ORDER BY read_sum DESC, hash ASC`. Der Hash ist inhaltsabhängig und damit unabhängig von Dateireihenfolge, Threads und Parquet-Lesereihenfolge. **vsearch-Vorbild (geprüft in `src/commands/sortbysize.cpp` und der Man-Page von `--sortbysize`/`--cluster_size`):** absteigende Abundanz, bei Gleichstand das Label in alphanumerischer Reihenfolge, zuletzt die Eingabereihenfolge. Unser Label ist der Hash, `ORDER BY read_sum DESC, hash ASC` entspricht also genau der vsearch-Regel (ohne den letzten Fallback, der bei eindeutigen Hashes nie greift). Prüfpunkt: DuckDB-VARCHAR-Vergleich ist byteweise, bei Hex-Hashes identisch zu ASCII-Ordnung.
Dazu: `sample_idx` per `ORDER BY sample`, `sequence_idx` per `ORDER BY hash` statt `row_number() OVER ()`.
- [x] Impl (zusätzlich Tie-Breaker bei `group_data.sequence_order` und im Log-Sort)
- [x] Test (groß, ohne Excel): Lauf 1 vs. Lauf 2 vs. Lauf mit `cores to use` = 2 → 0 Abweichungen in allen Tabellen; `sequence_order` entspricht exakt `(read_sum DESC, hash ASC)` für Sequenzen und Gruppen; klein: 2 Läufe (Kerne Standard / 1) identisch
- [x] Bench: kein Performance-Effekt messbar (siehe Protokoll)
- [x] Commit

**Neue Referenz ab hier:** `D:pscale_test_data\step_1_1\{large_1,small_1}`. Alle folgenden Schritte (außer 1.2/2.x mit gewollter Ergebnisänderung) müssen exakt diese Snapshots reproduzieren (`tests/compare_runs.py`). Gegenüber `main` geändert: 21 070 `sequence_order`, 35 Gruppenzuordnungen, 24 Repräsentanten (Folge der Tie-Auflösung).

### 1.2 Hash nach Normalisierung berechnen (Review 2.7)
`upper().strip()` **vor** `sha3_256`. Hinweis: ändert Hashes nur, wenn Eingaben Kleinbuchstaben/Whitespace enthalten.
- [x] Impl
- [x] Test: `tests/test_fasta_to_parquet.py` (`acgt` und `ACGT` → gleicher Hash; Header-Hash bleibt ohne Hashing unverändert); schlägt auf altem Code fehl, auf neuem grün
- [x] Regression groß: Snapshot identisch zu `step_1_1/large_1` (erwartet, da die Daten Header-Hashes nutzen)
- [x] Commit

---

## Phase 2 – Korrektheit

### 2.1 Gruppen-Parquet nutzt `sequence_data` (Review 2.1)
Reproduziert: Demo mit 4 Seqs / 3 Gruppen liefert 4 Zeilen, Zeile des Gruppenmitglieds ist NULL. Fix: überall `{sequence_type}_data` (Z. 371, 425), `sequence_ids` ebenfalls.
- [x] Impl (`sequence_ids`-Abfrage und finaler `FROM` nutzen `{sequence_type}_data`)
- [x] Test: `tests/check_read_tables.py <run_dir>` prüft Zeilenzahl, Zeilenreihenfolge == `sequence_order`, keine NULLs, Read-Summe je Sample. **Vorher** (echte Daten, groß): Gruppen-Parquet 27 964 statt 6 485 Zeilen, 21 479 Zeilen mit NULL. **Nachher:** alle 8 Prüfungen grün auf groß und klein; Snapshot unverändert gegenüber `step_1_1`
- [x] Commit

### 2.2 Zählertyp `HUGEINT` → `BIGINT` (Review 2.6)
Betrifft auch `group_data.read_sum` (im Demo bereits `210.0`). `::BIGINT` bei allen `SUM(...)`.
- [ ] Impl
- [ ] Test: `typeof` aller Zählspalten im Parquet und in den Tabellen ist `BIGINT`
- [ ] Commit

### 2.3 Chunking ohne Sortierlücke (Review 2.3)
Bereiche über `sequence_order BETWEEN a AND b` (arithmetisch, keine Python-Liste aller IDs). Wird mit 3.4 zusammen sauber, kann aber vorab einzeln gefixt werden, falls 3.4 später kommt.
- [ ] Impl
- [ ] Test: Anzahl Zeilen == Anzahl Sequenzen, keine Duplikate, bei absichtlich vertauschter Tabellenreihenfolge
- [ ] Commit

### 2.4 Log-Dateien (Review 2.2)
`file_nr` erhöhen, alte `Logfile_11_*` vorher löschen, Gruppen-**Hash** statt nur Indizes ausgeben. Schreiben mit DuckDB (siehe 3.1).
- [ ] Impl
- [ ] Test: Zeilen im Log == Zeilen in `group_mapping`, bei erzwungen kleinem Chunk mehrere Dateien
- [ ] Commit

### 2.5 Excel-Spaltenlimit (Review 2.5)
`chunksize = min(10_000_000 // sequence_count, 16_000)`; ≥ 1 sicherstellen.
- [ ] Impl
- [ ] Test: 20 000 Samples × 5 Seqs → mehrere Dateien, je ≤ 16 384 Spalten
- [ ] Commit

### 2.6 Samplenamen: Quoting und Kollisionen (Review 2.9)
Helper `quote_ident()` / `quote_literal()` für alle SQL-Stellen. Frühzeitiger, klarer Fehler bei Namen, die case-insensitiv kollidieren oder `hash`/`sequence` heißen.
- [ ] Impl
- [ ] Test: Samplenamen mit `'`, `"`, Leerzeichen, Punkten; `Sample1` vs. `sample1` → verständliche Fehlermeldung
- [ ] Commit

### 2.7 Alte Ausgaben aufräumen (Review 2.10)
Vor dem Schreiben `*_read_table_part_*.xlsx`, Parquet und `Logfile_11_*` des Projekts löschen.
- [ ] Impl
- [ ] Test: Lauf mit 3 Teilen, danach Lauf mit 1 Teil → nur 1 Teil übrig
- [ ] Commit

### 2.8 Eingaben und Fehlerbehandlung in `main` (Review 2.12, 2.13)
Keine `*.fasta.gz` → klare Meldung und Abbruch; ungültiger Schwellenwert → Exception bzw. Exit-Code ≠ 0; `mkdir(parents=True, exist_ok=True)`; Temp-Aufräumen robust (`shutil.rmtree(..., ignore_errors=True)`).
- [ ] Impl
- [ ] Test: leerer Ordner, Schwellenwert 0 / 1.5, vorhandener Temp-Ordner mit Resten
- [ ] Commit

### 2.9 vsearch-Aufruf (Review 2.11)
`check=True` bzw. Rückgabecode prüfen, klare Meldung wenn `vsearch` fehlt. Laut vsearch-Man-Page (Quelle geprüft): `--maxrejects` Default 32, nur bei `--maxaccepts 0` **und** `--maxrejects 0` wird die komplette DB durchsucht; Treffer-Kandidaten werden nach gemeinsamen Wörtern sortiert (Heuristik); `--qmask` Default `dust`. Unser Aufruf (`--maxaccepts 0`, maxrejects Default) kann also Treffer auslassen. **Testreihe**, ob das die Gruppen verändert:
`--maxrejects 0` vs. Default (32) und `--qmask/--dbmask none` vs. `dust`, auf mittlerem Datensatz und echten Daten.
- [ ] Impl (Fehlerprüfung)
- [ ] Test: vsearch nicht im PATH → saubere Fehlermeldung
- [ ] Experiment: Gruppen vergleichen (Anzahl Gruppen, Mapping-Differenz). Ergebnis hier notieren: ______
- [ ] Entscheidung (Parameter ändern ja/nein): ______
- [ ] Commit

---

## Phase 3 – Performance (nach Korrektheit, jeweils gegen Baseline messen)

### 3.1 Excel direkt über DuckDB (Wunsch: gesetzt)
`COPY (...) TO '*.xlsx' (FORMAT xlsx, HEADER true)` für Read-Tabellen und Log; kein pandas/xlsxwriter mehr.
**Offene Frage Memory:** ob der DuckDB-xlsx-Writer streamt oder die Sheet-Daten im RAM hält, weiß ich nicht. Der Benchmark entscheidet:
- 1 Mio. Zeilen × (2 + 10 Spalten) und 100 000 Zeilen × 5 000 Spalten, Peak-RSS messen, Vergleich mit altem Pfad
- Falls Peak-RSS hoch: Chunkgröße der Teil-Dateien runter, oder Fallback `xlsxwriter(constant_memory=True)` mit `fetchmany`-Streaming aus DuckDB
- Extension: `INSTALL excel`/`LOAD excel` mit klarer Meldung, wenn offline und nicht installiert
- [ ] Bench-Prototyp (vor Impl.)
- [ ] Impl
- [ ] Test: Werte Zelle für Zelle gleich wie Baseline (per openpyxl, `read_only=True`)
- [ ] Bench: Laufzeit + Peak-RSS vs. Baseline
- [ ] Commit

### 3.2 Kreuzprodukt entfernen (Review 3.2)
Pivot direkt auf der sparsen `*_read_count_data`, Nullen per `COALESCE(col, 0)` nach dem Pivot. `cross_product`, `read_table_data`, `filtered`, `pivot` als Tabellen entfallen. Excel-Teile filtern auf `sample_idx`-Bereiche.
Memory: nur DuckDB-intern (spillt), nichts in Python.
- [ ] Impl
- [ ] Test: Ausgabe identisch zu Baseline (Nullen statt NULL)
- [ ] Bench: Laufzeit, Peak-RSS, **Plattenverbrauch des Temp-Ordners**
- [ ] Commit

### 3.3 Gruppierung: O(N²)-Schleife ersetzen (Review 3.1)
**Memory-sicherer Entwurf** (nicht „alles im RAM“):
- Match-Tabelle bleibt in DuckDB. Abfrage: `... WHERE target_order > query_order ORDER BY query_order, target_order` (DuckDB sortiert out-of-core), per `fetchmany` gestreamt.
- In Python nur ein `numpy`-int32-Array `group_of[order]` (Größe = Anzahl Sequenzen, ~4 MB pro 1 Mio.).
- Ein Durchlauf: Für jede Query in Reihenfolge: ist sie noch nicht zugeordnet → sie ist Repräsentant; alle Targets mit noch freiem Slot bekommen diese Gruppe. Am Ende: jede nicht zugeordnete Sequenz ist eine eigene Gruppe.
- Semantik identisch zum bisherigen Code (Treffer auf Sequenzen mit kleinerer Order werden ignoriert, weil die schon entschieden sind).
- [ ] Alte Implementierung als Test-Orakel in `tests/` sichern (nach 1.1, damit beide dieselbe Order nutzen)
- [ ] Impl
- [ ] Test: neuer Algorithmus == altes Orakel auf Zufallsdaten (viele Seeds, mit Ties, ohne Treffer, erste Sequenz ohne Treffer, asymmetrische Treffer), identisches `group_mapping`
- [ ] Test: Lauf-zu-Lauf-Stabilität (Phase 1.1)
- [ ] Bench: Laufzeit und Peak-RSS bei 100k / 1 Mio. Sequenzen
- [ ] Commit

### 3.4 Parquet nach `sequence_order`-Bereichen schreiben (Review 3.4)
Jeder Chunk joint Hash/Sequenz selbst und schreibt in Reihenfolge; kein globaler `ORDER BY` über breite Zeilen mehr. Dateien werden danach in Reihenfolge zu einer Ausgabedatei zusammengeführt (Reihenfolgetest!).
- [ ] Impl
- [ ] Test: identisch zu Baseline inkl. Zeilenreihenfolge
- [ ] Bench: Laufzeit, Peak-RSS
- [ ] Commit

### 3.5 `build_read_store` schlanker (Review 3.5)
`sample_idx` schon in Python vergeben und in die Parquet-Dateien schreiben (kein String-Join auf `sample`). `sequence_idx` per Join auf `hash`; vorher Eindeutigkeit `hash → sequence` prüfen (Fehler, wenn nicht eindeutig). `ALTER` + `UPDATE` durch ein einziges `CREATE TABLE AS` mit Aggregation ersetzen.
- [ ] Impl
- [ ] Test: Tabellen identisch zu Baseline (bis auf Index-Vergabe, dort gegen 1.1 prüfen)
- [ ] Bench
- [ ] Commit

### 3.6 `fasta_to_parquet` ohne pandas (Review 3.6)
Listen → `pyarrow` direkt (falls pyarrow ohnehin Abhängigkeit von pandas/DuckDB ist; sonst DuckDB `from_arrow`/`executemany` prüfen).
- [ ] Impl
- [ ] Test: Parquet-Inhalt identisch
- [ ] Bench
- [ ] Commit

### 3.7 DuckDB-Konfiguration (Review 7)
`threads` aus `cores to use`, `memory_limit` (Anteil des freien RAMs), `temp_directory` auf Projekt-Temp, `preserve_insertion_order=false` für große `COPY`.
- [ ] Impl
- [ ] Bench: große Daten mit künstlich kleinem `memory_limit` laufen durch (Spill-Test)
- [ ] Commit

---

## Phase 4 – Kleinigkeiten (ein Sammel-Commit)
- [ ] Stabile Spaltenreihenfolge statt Tuple-Unpacking aus `SELECT *`
- [ ] Ungenutzte Variablen entfernen (`sequence_counter`, Zuweisung bei `COPY`, doppelte `pbar.update`)
- [ ] Docstrings korrigieren (`to_excel` entfernen), veraltetes „hdf store“-Kommentar
- [ ] Teilenummern vereinheitlichen (Excel und Log beginnen bei 1)
- [ ] Typ-Hints `str` → `Path`
- [ ] Logging: einheitliche Hilfsfunktion statt wiederholtem `datetime.now().strftime(...)`

---

## Abschluss
- [ ] Voller Lauf auf realem Datensatz: Ausgaben gegen Baseline vergleichen (bekannte, gewollte Unterschiede: Gruppen-Parquet, Log, BIGINT, Tie-Reihenfolge)
- [ ] `tests/bench_results.md`: Vorher/Nachher-Tabelle
- [ ] Manual/README anpassen, falls Ausgabedateien oder Verhalten sich ändern
- [ ] PR gegen `main`

## Benchmark-Protokoll (wird je Schritt ergänzt)

Alle Läufe mit den Standard-Settings aus `create_project` (`cores to use` = CPU-Kerne − 2, Schwelle 0.97).

| Schritt | Datensatz | Zeit alt | Zeit neu | Peak-RSS alt | Peak-RSS neu | Notiz |
|---|---|---|---|---|---|---|
| 0 Baseline `d6bb6fa` | small | 7,4 s / 6,9 s | – | 2,66 GB / 2,64 GB | – | 2 Läufe, Ergebnisse identisch |
| 0 Baseline `d6bb6fa` | large | 303 s / 307 s | – | 7,0 GB / 7,0 GB | – | 2 Läufe, Ergebnisse **nicht** identisch (s. u.) |

| 1.1 Tie-Breaker | large, `--skip-excel` | 303 s (mit Excel) | 156 s / 151 s / 151 s | 7,0 GB | 6,8 / 6,8 / 4,6 GB | Zeiten nicht vergleichbar (Excel aus); 3. Lauf mit 2 Kernen, gleiches Ergebnis |
| 1.1 Tie-Breaker | small, `--skip-excel` | 7 s | 4,2 s / 4,8 s | 2,6 GB | – | 2. Lauf mit 1 Kern |

**Befund Baseline large (zwei Läufe von unverändertem `main`):** 21 006 von 27 964 Sequenzen haben ein anderes `sequence_order`, 20 Sequenzen landen in einer anderen Gruppe, 12 Gruppenrepräsentanten unterscheiden sich, 15 Zeilen in `group_read_counts` weichen ab. `read_counts` identisch. Nur 6 965 verschiedene Read-Summen bei 27 964 Sequenzen (viele Ties). → Tie-Bug (Schritt 1.1) ist auf echten Daten nachgewiesen.
