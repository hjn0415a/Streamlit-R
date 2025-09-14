import shutil
from pathlib import Path

import streamlit as st

from src.common.common import reset_directory

def save_uploaded_excel(uploaded_files: list[bytes]) -> None:
    """
    Saves uploaded Excel files to the excel directory.
    """
    excel_dir = Path(st.session_state.workspace, "Excel-files")

    if st.session_state.location == "online":
        uploaded_files = [uploaded_files]

    if not uploaded_files:
        st.warning("Upload some Excel/CSV files first.")
        return

    for f in uploaded_files:
        if f.name not in [f.name for f in excel_dir.iterdir()] and f.name.endswith((".xlsx", ".csv")):
            with open(Path(excel_dir, f.name), "wb") as fh:
                fh.write(f.getbuffer())
    st.success("Successfully added uploaded Excel/CSV files!")

def copy_local_excel_files_from_directory(local_excel_directory: str, make_copy: bool = True) -> None:
    """
    Copies local Excel files from a specified directory to the Excel directory.
    """
    excel_dir = Path(st.session_state.workspace, "excel-files")
    valid_exts = (".xlsx", ".csv")

    files = [f for f in Path(local_excel_directory).iterdir() if f.suffix.lower() in valid_exts]

    if not files:
        st.warning("No Excel/CSV files found in specified folder.")
        return

    for f in files:
        if make_copy:
            shutil.copy(f, Path(excel_dir, f.name))
        else:
            external_files = Path(excel_dir, "external_files.txt")
            if not external_files.exists():
                external_files.touch()
            with open(external_files, "a") as f_handle:
                f_handle.write(f"{f}\n")

    st.success("Successfully added local Excel/CSV files!")

def remove_selected_excel_files(to_remove: list[str], params: dict) -> dict:
    """
    Removes selected Excel files from the excel directory.
    """
    excel_dir = Path(st.session_state.workspace, "excel-files")

    for f in to_remove:
        Path(excel_dir, f).unlink(missing_ok=True)

    for k, v in params.items():
        if isinstance(v, list) and any(f in v for f in to_remove):
            params[k] = [item for item in v if item not in to_remove]

    st.success("Selected Excel/CSV files removed!")
    return params

def remove_all_excel_files(params: dict) -> dict:
    """
    Removes all Excel files from the excel directory.
    """
    from src.common.common import reset_directory

    excel_dir = Path(st.session_state.workspace, "excel-files")
    reset_directory(excel_dir)

    for k, v in params.items():
        if "xlsx" in k and isinstance(v, list):
            params[k] = []

    st.success("All Excel/CSV files removed!")
    return params