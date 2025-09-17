import os
import streamlit as st
import pandas as pd

from src.common.common import page_setup

params = page_setup()

data_dir = "/data"

st.markdown("## Upload CSV File")

uploaded_file = st.file_uploader(
        "Upload a CSV file",
        type=["csv"],
        accept_multiple_files=False
)

submit = st.button("Submit")

if submit:
    if uploaded_file is not None:
        save_path = os.path.join(data_dir, uploaded_file.name)
        with open(save_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        st.success(f"File saved to {save_path}")

        try:
            df = pd.read_csv(save_path)
            st.dataframe(df)
        except Exception as e:
            st.error(f"Error reading CSV: {e}")
    else:
        st.error("Please upload a CSV file first.")
