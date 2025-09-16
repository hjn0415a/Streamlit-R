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
        "FullseePathway": [
             st.Page(Path("content", "rcode", "quickstart.py"), title="User Guide", icon="📃"),
             st.Page(Path("content", "rcode", "volcano_plot.py"), title="File Upload", icon="📁"),
        ],
        "Basic analysis": [
             st.Page(Path("content", "rcode", "temp.py"), title="Heatmap", icon="🌡️"),
             st.Page(Path("content", "rcode", "temp.py"), title="Volcano, Enhanced", icon="🌋"),
             st.Page(Path("content", "rcode", "temp.py"), title="PCA", icon="📉"),
        ],
        "GO Pathway analysis": [
             st.Page(Path("content", "rcode", "temp.py"), title="Enrichplot", icon="📈"),
             st.Page(Path("content", "rcode", "temp.py"), title="Cnetplot", icon="🕸️"),
             st.Page(Path("content", "rcode", "temp.py"), title="Emapplot", icon="📊"),
             st.Page(Path("content", "rcode", "temp.py"), title="Ridgeplot", icon="🕸️"),
             st.Page(Path("content", "rcode", "temp.py"), title="Heatmaplike functional classification", icon="📋"),
        ],
        "Gene-gene Interaction": [
             st.Page(Path("content", "rcode", "temp.py"), title="Enrichplot", icon="📈"),
        ],
        "Enrichkegg": [
             st.Page(Path("content", "rcode", "temp.py"), title="Pathview", icon="🔎"),
        ]
    }

    pg = st.navigation(pages)
    pg.run()