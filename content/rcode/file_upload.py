import os
import streamlit as st
import pandas as pd

data_dir = "/data"

st.markdown("## Upload CSV File")

# Form 생성
with st.form(key="upload_form"):
    uploaded_file = st.file_uploader(
        "Upload a CSV file",
        type=["csv"],
        accept_multiple_files=False
    )
    submit_button = st.form_submit_button(label="Submit")

# 폼 제출 시 처리
if submit_button and uploaded_file is not None:
    # 파일 저장
    save_path = os.path.join(data_dir, uploaded_file.name)
    with open(save_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    st.success(f"File saved to {save_path}")

    try:
        # CSV 읽기
        df = pd.read_csv(save_path)
        st.markdown("### Uploaded Table")
        st.dataframe(df)
    except Exception as e:
        st.error(f"Error reading CSV: {str(e)}")
