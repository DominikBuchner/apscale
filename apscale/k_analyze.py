import streamlit as st
import apscale, subprocess, sys
from pathlib import Path


def main(project=Path.cwd()):
    # find installation of apscale
    st.title("Apscale analysis module")
    project = Path(sys.argv[1])
    st.write(project)


if __name__ == "__main__":
    main()
