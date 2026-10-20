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
- [x] Impl (`::BIGINT` an den 5 `SUM`-Stellen: `read_sum`, beide `PIVOT`s, `group_data`, `group_read_count_data`)
- [x] Test: `check_read_tables.py` prüft Spaltentypen. **Vorher:** Parquet-Zählspalten `DOUBLE`, `group_read_count_data.read_count` `DOUBLE`/`HUGEINT`. **Nachher:** alle `BIGINT`, 12 Prüfungen grün auf groß; Snapshot-Inhalt unverändert gegenüber `step_1_1`
- [x] Commit

### 2.3 Chunking ohne Sortierlücke (Review 2.3)
Bereiche über `sequence_order BETWEEN a AND b` (arithmetisch, keine Python-Liste aller IDs). Wird mit 3.4 zusammen sauber, kann aber vorab einzeln gefixt werden, falls 3.4 später kommt.
- [x] Impl: Chunk-Bereiche per SQL (`row_number() OVER (ORDER BY sequence_idx) // Chunkgröße` → min/max je Chunk), keine Python-Liste mehr; Chunkgrößen als Modulkonstanten `MAX_PIVOT_CELLS`/`MIN_ROWS_PER_CHUNK` (Verhalten unverändert, testbar)
- [x] Test: `tests/test_chunking.py` (reale `small_dataset`, erzwungen 7 Zeilen pro Chunk, Sequenzen/Gruppen × Tabelle unvertauscht/vertauscht). **Vorher:** 3 von 4 Fällen rot, u. a. 1 503 statt 120 Zeilen bei den Gruppen **ohne** Vertauschen (`group_data` ist durch `GROUP BY` nicht nach `sequence_idx` sortiert), 22 227 statt 456 bei vertauschten Sequenzen. **Nachher:** 4/4 grün; groß: 12 Prüfungen grün, Snapshot unverändert
- [x] Commit

### 2.4 Log-Dateien (Review 2.2)
`file_nr` erhöhen, alte `Logfile_11_*` vorher löschen, Gruppen-**Hash** statt nur Indizes ausgeben. Schreiben mit DuckDB (siehe 3.1).
- [x] Impl: `fetchmany(LOG_ROWS_PER_FILE)` (200 000), `file_nr` wird erhöht, alte `Logfile_11_*` werden vorher gelöscht, Spalten jetzt `group_hash`, `sequence_hash`, `read_sum` (interne Indizes entfernt). Schreiben vorerst noch über pandas/xlsxwriter (Umstellung in 3.1)
- [x] Test: `tests/test_write_log.py` (mehrere Dateien bei kleinem Limit, Zeilenzahl == `group_mapping`, Zuordnung stimmt mit DB überein, Sortierung, Altdatei wird entfernt). **Vorher** (Demo, 600-fach vervielfachte Zuordnung, 273 600 Zeilen): nur 1 Datei mit 68 800 Zeilen, 204 800 Zeilen verloren. **Nachher:** 200 000 + 73 600 in 2 Dateien. Vollständiger Lauf groß mit Excel: Log 27 964 Zeilen, Snapshot unverändert
- [x] Commit

### 2.5 Excel-Spaltenlimit (Review 2.5)
`chunksize = min(10_000_000 // sequence_count, 16_000)`; ≥ 1 sicherstellen.
- [x] Impl: `chunksize = max(1, min(EXCEL_MAX_CELLS // sequence_count, EXCEL_MAX_COLUMNS - 2))` (Konstanten am Modulanfang). Zeilenlimit-Schwelle (1 000 000) bewusst unverändert
- [x] Test: `tests/test_excel_limits.py` (synthetisch, 17 000 Samples × 3 Sequenzen, dort ohne echte Daten möglich). **Vorher:** `ValueError: This sheet is too large! Your sheet size is: 3, 17002 Max sheet size is: 1048576, 16384` – der ganze Lauf bricht ab, auch das Parquet entsteht nicht mehr. **Nachher:** ≥ 2 Teile, je ≤ 16 384 Spalten, alle Samples genau einmal vorhanden. Groß: Snapshot und 12 Prüfungen unverändert
- [x] Commit

### 2.6 Samplenamen: Quoting und Kollisionen (Review 2.9)
Helper `quote_ident()` / `quote_literal()` für alle SQL-Stellen. Frühzeitiger, klarer Fehler bei Namen, die case-insensitiv kollidieren oder `hash`/`sequence` heißen.
- [x] Impl: `quote_literal()`/`quote_ident()` für alle Samplenamen **und alle Pfade** in SQL (`ATTACH`, `read_parquet`, `read_csv`, `COPY ... TO`); `check_sample_names()` (reservierte Namen `hash`/`sequence`/`sequence_idx`, Namen die sich nur in Groß-/Kleinschreibung unterscheiden) wird früh in `main` (aus den Dateinamen) und in `generate_read_table` aufgerufen
- [x] Test: `tests/test_sample_names.py` (Namen mit `'`, `"`, Leerzeichen, Punkt, `;`, Umlauten, `[brackets]*`; Projektordner mit Apostroph; Datei `O'Brien.fasta.gz` durch `fasta_to_parquet`; Kollisionen). **Vorher:** `Parser Error: syntax error at or near "s"` schon bei Projektordnern mit Apostroph (Pfade im SQL), Kollisionsnamen liefen ohne Fehlermeldung durch. **Nachher:** 7/7 grün, gesamte Suite 15/15; groß und klein: Snapshot unverändert
- [x] Commit

### 2.7 Alte Ausgaben aufräumen (Review 2.10)
Vor dem Schreiben `*_read_table_part_*.xlsx`, Parquet und `Logfile_11_*` des Projekts löschen.
- [x] Impl: `remove_old_read_tables(folder, project_name, sequence_type=None)` löscht `*_read_table_part_*.xlsx` und `*_read_table.parquet.snappy` des Projekts; Aufruf am Anfang von `generate_read_table` (für den jeweiligen Typ) und in `main` für beide Typen (auch wenn „generate read table“ ausgeschaltet ist). Das Log wird seit 2.4 in `write_log` bereinigt
- [x] Test: `tests/test_stale_outputs.py` (4 Teile → danach 1 Teil: **vorher** blieben `part_1..3` liegen, **nachher** nur `part_0`; fremde Dateien bleiben unberührt; Parquet und beide Typen). End-to-end: `main` auf `small_dataset` mit abgelegten Alt-Dateien (Excel-Teil, Parquet, Log) → alle weg bzw. überschrieben. Suite 18/18, groß: Snapshot unverändert
- [x] Commit

### 2.8 Eingaben und Fehlerbehandlung in `main` (Review 2.12, 2.13)
Keine `*.fasta.gz` → klare Meldung und Abbruch; ungültiger Schwellenwert → Exception bzw. Exit-Code ≠ 0; `mkdir(parents=True, exist_ok=True)`; Temp-Aufräumen robust (`shutil.rmtree(..., ignore_errors=True)`).
- [x] Impl: `ValueError` bei Schwellenwert ∉ (0, 1] (statt `print` + `sys.exit()`), `ValueError` wenn kein vorheriger Schritt aktiv ist (`choose_input` liefert `None`), `FileNotFoundError` ohne `*.fasta.gz`; Ordner mit `mkdir(parents=True, exist_ok=True)`; Temp-Ordner wird zu Beginn komplett neu angelegt und am Ende mit `shutil.rmtree` entfernt; ungenutzte Imports `os`/`sys` entfernt
- [x] Test: `tests/test_main_inputs.py` (6 Fälle: Schwellenwerte 0 / 1.5 / −0.2, kein Input, Reste im Temp-Ordner, fehlender `11_read_table`-Ordner). **Vorher:** alle rot – `SystemExit` ohne Fehlercode, kryptischer DuckDB-IOError bei fehlendem Input, **Reste im Temp-Ordner (`sequence_read_count_data.duckdb`) brechen den ganzen Lauf**, `FileNotFoundError` bei fehlendem Ordner. **Nachher:** alle grün, Suite 24/24; groß und klein: Snapshot unverändert
- [x] Commit

### 2.9 vsearch-Aufruf (Review 2.11)
`check=True` bzw. Rückgabecode prüfen, klare Meldung wenn `vsearch` fehlt. Laut vsearch-Man-Page (Quelle geprüft): `--maxrejects` Default 32, nur bei `--maxaccepts 0` **und** `--maxrejects 0` wird die komplette DB durchsucht; Treffer-Kandidaten werden nach gemeinsamen Wörtern sortiert (Heuristik); `--qmask` Default `dust`. Unser Aufruf (`--maxaccepts 0`, maxrejects Default) kann also Treffer auslassen. **Testreihe**, ob das die Gruppen verändert:
`--maxrejects 0` vs. Default (32) und `--qmask/--dbmask none` vs. `dust`, auf mittlerem Datensatz und echten Daten.
- [x] Impl (Fehlerprüfung: `RuntimeError` bei fehlendem vsearch bzw. Rückgabecode ≠ 0 mit stderr; `VSEARCH_EXTRA_ARGS` für Experimente)
- [x] Test: vsearch nicht im PATH / vsearch schlägt fehl → klare Fehlermeldung (`tests/test_main_inputs.py`)
- [x] Experiment (`tests/experiments/vsearch_params.py`, groß, Schwelle 0.97): Maskierung (`--qmask/--dbmask none`) ändert nichts (7/6 Paare, 0 Gruppen). `--maxrejects 0` findet 65 727 zusätzliche Paare (+6,8 %, alle bei 97,0–98,5 % Identität), 6 470 statt 6 485 Gruppen, 441 Sequenzen in anderer Gruppe; vsearch-Zeit ca. 11 min statt < 1 min
- [x] Entscheidung: **`--maxrejects` bleibt beim Standard 32**, Maskierung bleibt Standard. Siehe aber 3.3: vsearch selbst (`--cluster_size`) vermeidet das Problem
- [x] Commit

### 3.1 Excel direkt über DuckDB (Wunsch: gesetzt)
`COPY (...) TO '*.xlsx' (FORMAT xlsx, HEADER true)` für Read-Tabellen und Log; kein pandas/xlsxwriter mehr.
**Offene Frage Memory:** ob der DuckDB-xlsx-Writer streamt oder die Sheet-Daten im RAM hält, weiß ich nicht. Der Benchmark entscheidet:
- 1 Mio. Zeilen × (2 + 10 Spalten) und 100 000 Zeilen × 5 000 Spalten, Peak-RSS messen, Vergleich mit altem Pfad
- Falls Peak-RSS hoch: Chunkgröße der Teil-Dateien runter, oder Fallback `xlsxwriter(constant_memory=True)` mit `fetchmany`-Streaming aus DuckDB
- Extension: `INSTALL excel`/`LOAD excel` mit klarer Meldung, wenn offline und nicht installiert
- [x] Bench-Prototyp (`tests/experiments/excel_writer_bench.py`, größte reale Dateigröße ≈ 10 Mio. Zellen, jede Messung in eigenem Prozess):

| Fall | pandas + xlsxwriter (alt) | DuckDB `COPY ... (FORMAT xlsx)` | xlsxwriter `constant_memory` + `fetchmany` |
|---|---|---|---|
| lang: 800 000 Zeilen × 12 Spalten | 88 s / **2 203 MB** | **13,8 s / 45 MB** | 49 s / 91 MB |
| breit: 10 000 Zeilen × 1 002 Spalten | 76 s / **1 688 MB** | **7,9 s / 135 MB** | 34 s / 391 MB |

  → der DuckDB-Writer **streamt** (Speicher unabhängig von der Dateigröße) und ist 6–10× schneller. Inhalt geprüft: Kopfzeile, Zeilenzahl, erste Zeile, Zelltypen (`str`/`int`) und Sheet-Name (`Sheet1`) identisch zur pandas-Variante. Entscheidung: DuckDB-Writer, kein Fallback nötig.
- [x] Impl: Hilfsfunktion `write_excel(connection, query, file)` (lädt die `excel`-Extension, installiert sie bei Bedarf, klare Fehlermeldung wenn offline); Read-Tabellen und Log schreiben per `COPY`, kein `.df()`/pandas/xlsxwriter mehr in diesem Modul. Das Log wird einmal nummeriert (`TEMP TABLE log_rows`) und in Zeilenbereichen geschrieben
- [x] Test: `tests/test_excel_content.py` (Sequenzen und Gruppen, viele Teile: Hash, Sequenz, Reihenfolge und Zählwerte == Parquet, jedes Sample genau einmal), bestehende Excel-/Log-Tests weiter grün, Suite 33/33. Echte Daten `large`: Excel-Teile == Parquet (Sequenzen: 2 Teile, 27 964 Zeilen; Gruppen: 1 Teil, 6 470 Zeilen), Snapshot unverändert gegenüber `step_3_3`
- [x] Bench (`large`, **mit** Excel, gesamte Pipeline): **43 s statt 303 s**; Peak-RSS unverändert 7,0 GB (also nicht vom Excel-Schreiben, offen für 3.2/3.7)
- [x] Commit

### 3.2 Kreuzprodukt entfernen (Review 3.2)
Pivot direkt auf der sparsen `*_read_count_data`, Nullen per `COALESCE(col, 0)` nach dem Pivot. `cross_product`, `read_table_data`, `filtered`, `pivot` als Tabellen entfallen. Excel-Teile filtern auf `sample_idx`-Bereiche.
Memory: nur DuckDB-intern (spillt), nichts in Python.
- [x] Impl: `pivot_query()` pivotiert direkt die sparse `*_read_count_data` (Excel-Teile: pro Sample-Chunk als Unterabfrage ohne Zwischentabellen; Parquet: pro `sequence_idx`-Bereich), `zero_filled_columns()` setzt fehlende Werte per `COALESCE(..., 0)::BIGINT`. Entfernt: `cross_product`, `read_table_data`, `filtered`
- [x] Test: Suite 33/33; `large`: beide Parquet-Read-Tabellen **zellgleich** zu Schritt 3.1 (0 Zeilen Unterschied in beide Richtungen), Excel == Parquet weiterhin getestet, Snapshot unverändert, 12 Prüfungen grün
- [x] Bench (`large`, mit Excel): **24 s statt 43 s**, Peak-RSS **3,2 GB statt 7,0 GB**, Hauptprozess maximal 744 MB (vorher 4,5 GB). Speicherprofil nach Phasen mit `tests/phase_memory.py`
- [x] Befund aus dem Profil: **die restlichen ca. 2,4 GB sind 23 joblib-Workerprozesse** (loky), die nach `fasta_to_parquet` weiterleben (siehe 3.6b)
- [x] Commit

### 3.3 Gruppierung: eigene O(N²)-Schleife ersetzen (Review 3.1)
**Neuer Befund aus dem Vergleich mit dem vsearch-Quellcode (`src/core/cluster.cpp`) und einem Test auf den echten Daten:**
- vsearch clustert **online gegen die bisherigen Zentroide** (Datenbank = nur Zentroide), unsere Methode sucht **alle gegen alle** und fügt danach in Python zusammen. Bei `maxrejects 32` verbraucht der All-vs-all-Aufruf sein Budget an Nicht-Zentroiden und verpasst Nachbarn nahe der Schwelle (Ursache plausibel, nicht separat nachgewiesen).
- vsearch sortiert Eingaben intern nach `(size desc, Label alphanumerisch, Eingabereihenfolge)`. Unsere Reihenfolge `(read_sum desc, hash asc)` aus 1.1 ist damit identisch.
- `vsearch --cluster_size in.fasta --id 0.97 --sizeorder --maxaccepts 0 --uc out.uc` (Header `hash;size=read_sum`, **maxrejects Standard 32**) liefert auf `large` **exakt** dieselbe Zuordnung wie unser exhaustiver All-vs-all-Lauf (6 470 Gruppen, 0 Abweichungen), in **2 s** (1 Thread: 19 s) statt rund 110 s. Ergebnis ist unabhängig von der Threadzahl (1 / 3 / alle: identische `.uc`). Ohne `--sizeorder` wählt vsearch den nächsten statt des häufigsten Zentroids → andere Semantik (2 107 Abweichungen), nicht verwenden.
- Gegenüber dem heutigen Stand ändert das 441 Zuordnungen und 15 Gruppen (die heute verpassten Treffer). **Das ist eine Ergebnisänderung und braucht eine Freigabe.**

Vorschlag (statt Python-Schleife): `generate_sequence_groups` schreibt `hash;size=read_sum` als FASTA, ruft `vsearch --cluster_size ... --uc` auf, liest die `.uc` (Typ `S` = Zentroid, `H` = Mitglied) per DuckDB und baut `group_mapping` daraus. Entfällt: Matchfile, `group_lookup_table`, `groups_found`, Python-Schleife, tqdm. Memory: O(Sequenzen) in vsearch, nichts in Python.
- [x] Freigabe der Ergebnisänderung (Nutzer: „Mir scheint es sinnvoll die vsearch-Variante nachzubauen“; Mapping Sequenz → OTU bleibt erhalten)
- [x] Alte Implementierung als Test-Orakel gesichert: `tests/oracle_grouping.py` (All-vs-all, `--maxrejects 0`)
- [x] Impl: `generate_sequence_groups` schreibt `hash;size=read_sum` als FASTA, ruft **einmal** `vsearch --cluster_size --sizeorder --maxaccepts 0 --minseqlength 1 --uc` auf (maxrejects Standard 32) und baut `group_mapping` aus der `.uc` (`S` = Zentroid, `H` = Mitglied). Prüfungen: Hashes eindeutig, jede Sequenz im Mapping, leere Eingabe. `--minseqlength 1`, damit vsearch keine kurzen Sequenzen stillschweigend verwirft. Entfernt: Matchfile, `group_lookup_table`, `groups_found`, Python-Schleife, `tqdm`
- [x] Test: `tests/test_grouping.py` (neues Mapping == Orakel bei Schwelle 0.97 und 0.9, Threadzahl 1 vs. 3 identisch, Konsistenz der Gruppentabellen, Einzelsequenz, leere Eingabe), Suite 31/31. **Echte Daten `large`:** neues Mapping == erschöpfender Referenzlauf des alten Verfahrens (**0 Abweichungen**, 6 470 Gruppen); gegenüber `step_1_1` wie angekündigt 441 Sequenzen / −15 Gruppen
- [x] Bench (`large`, ohne Excel, gesamte Pipeline): **ca. 30 s statt 138–156 s**, Gruppierung selbst ca. 2 s statt ca. 110 s. Peak-RSS unverändert 7,0 GB (kommt aus anderen Schritten, siehe 3.1/3.2)
- [x] Neue Referenz-Snapshots: `D:pscale_test_data\step_3_3\{large,small}`; ab hier gegen diese vergleichen
- [x] Commit

**Nachtrag Memory (vsearch selbst, `tests/experiments/vsearch_memory.py`, synthetisch skalierte Daten):** Peak-RSS des vsearch-Prozesses bei 28k / 100k Sequenzen: neu (`cluster_size`) 54 / 154 MB, alt (`usearch_global` alle gegen alle) 61 / 156 MB, also etwa linear und kein Vorteil für die alte Variante; die alte DuckDB-Schleife war nur beim Zusammenführen speichersicher, nicht beim vsearch-Lauf. Der 300k-Lauf wurde abgebrochen (Frage beantwortet). Offen/optional: `--maxaccepts 100` statt `0` war bei 100k gut dreimal schneller (8,8 s statt 28,9 s); Gleichheit der Ergebnisse noch **nicht** geprüft, daher bisher bei `0`.

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

### 3.6b joblib-Worker nach dem Einlesen beenden (neu, aus dem Speicherprofil)
Die loky-Workerprozesse bleiben nach `Parallel(...)` bestehen (wiederverwendbarer Executor) und belegen für den Rest des Laufs ca. 100 MB je Prozess (bei dir 23 Prozesse ≈ 2,4 GB), obwohl sie nur beim Einlesen der FASTA-Dateien gebraucht werden.
- [ ] Impl: Executor nach dem Einlesen herunterfahren (`get_reusable_executor().shutdown(wait=True)`) oder Parallel mit eigenem Backend-Kontext
- [ ] Test/Bench: Peak-RSS `large` vor/nach, Ergebnis unverändert
- [ ] Commit

### 3.7 DuckDB-Konfiguration (Review 7)
`threads` aus `cores to use`, `memory_limit` (Anteil des freien RAMs), `temp_directory` auf Projekt-Temp, `preserve_insertion_order=false` für große `COPY`.
- [ ] Impl
- [ ] Bench: große Daten mit künstlich kleinem `memory_limit` laufen durch (Spill-Test)
- [ ] Commit

### 3.8 DuckDB-Version (ganz am Ende der Performance-Phase)
Aktuell: `setup.py` verlangt `duckdb>=1.3.1`, getestet wird nur mit 1.5.2 (PyPI: 1.5.6). Bewusst **nicht mitten in den Umbau**, damit Benchmarks nicht durch Versionswechsel verfälscht werden.
- [ ] Gesamte Suite + Harness auf groß/klein mit 1.5.2 und der dann aktuellen 1.5.x laufen lassen (Laufzeit, Peak-RSS, Snapshots identisch)
- [ ] Prüfen, welche Mindestversion wir wirklich brauchen (xlsx-`COPY` aus 3.1, `PIVOT`, `fetchmany`) und `setup.py` auf die niedrigste **getestete** Version setzen
- [ ] Dateiformat der `*_read_data_storage.duckdb` wird von `Apscale_analyze.py` weiterverwendet: prüfen, dass eine mit neuerer Version erzeugte Datei mit der Mindestversion lesbar bleibt, sonst Mindestversion anheben
- [ ] Commit
- Nicht geplant: GPU-Erweiterung `gpudb` (Apple Metal / Linux CUDA, kein Windows, siehe Antwort im Chat)

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
| 2.4 Log | large, mit Excel | 303 s / 307 s | 319 s | 7,0 GB | 7,8 GB | volle Pipeline mit Excel; Abweichung im Rahmen der Messschwankung bzw. Rest-Overhead, nicht untersucht |
| 3.1 Excel per DuckDB | large, **mit** Excel | 303 s / 307 s | 43 s | 7,0 GB | 7,0 GB | Excel-Schreiben selbst: 6–10× schneller, Speicher −90 % im Prototyp; Gesamtspeicher kommt woanders her |
| 3.2 Kreuzprodukt entfernen | large, mit Excel | 43 s | 24 s | 7,0 GB | 3,2 GB | Hauptprozess 4,5 GB → 0,74 GB; Rest sind joblib-Worker (3.6b) |
| 3.3 Gruppierung per `vsearch --cluster_size` | large, `--skip-excel` | 138–156 s | 30 s | 6,8–7,0 GB | 7,0 GB | Gruppierung 110 s → 2 s; Mapping == exhaustiver Referenzlauf |
| 1.1 Tie-Breaker | small, `--skip-excel` | 7 s | 4,2 s / 4,8 s | 2,6 GB | – | 2. Lauf mit 1 Kern |

**Befund Baseline large (zwei Läufe von unverändertem `main`):** 21 006 von 27 964 Sequenzen haben ein anderes `sequence_order`, 20 Sequenzen landen in einer anderen Gruppe, 12 Gruppenrepräsentanten unterscheiden sich, 15 Zeilen in `group_read_counts` weichen ab. `read_counts` identisch. Nur 6 965 verschiedene Read-Summen bei 27 964 Sequenzen (viele Ties). → Tie-Bug (Schritt 1.1) ist auf echten Daten nachgewiesen.
