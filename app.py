import streamlit as st
from pathlib import Path
import json
# For some reason the windows version only works if this is imported here
import pyopenms

if "settings" not in st.session_state:
        with open("settings.json", "r") as f:
            st.session_state.settings = json.load(f)

if __name__ == '__main__':
    pages = {
        "QuantMS": [
             st.Page(Path("content", "quantms", "quickstart.py"), title="Quickstart", icon="👋"),
             st.Page(Path("content", "quantms", "volcano_plot.py"), title="R-Integration", icon="📃"),
        ]
    }

    pg = st.navigation(pages)
    pg.run()