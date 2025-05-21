import streamlit as st
import apscale, subprocess, sys
from pathlib import Path
import pandas as pd


def main(project=Path.cwd()):
    # configure the sidebar
    st.title("Apscale analysis module")
    try:
        project = Path(sys.argv[1])
    except IndexError:
        project = project

    st.write(f"Current project: **{project.name}**")


if __name__ == "__main__":
    main()
