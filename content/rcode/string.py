import os
import subprocess
import tempfile
import shutil
import streamlit as st

from src.common.common import page_setup
params = page_setup()

st.title("STRING Network Analysis Dashboard")

# ----------------- Main Tabs -----------------
main_tabs = st.tabs(["📊 STRING Network"])
with main_tabs[0]:
    # ----------------- Sub Tabs -----------------
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Configure -----------------
    with configure_tab:
        # 고정값으로 정의
        input_root = "/data/Deg"
        output_dir = "/data/STRING"
        os.makedirs(output_dir, exist_ok=True)

        combo_file = "combo_names.rds"
        taxon_id   = st.selectbox("Taxon ID", [9606, 10090], index=0, help="9606=Human, 10090=Mouse")
        cutoff     = st.slider("Confidence cutoff", 0.0, 1.0, 0.5, 0.05)
        limit      = st.number_input("Max interactions per gene (0=all)", value=0, step=1)
    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run STRING network generation"):
            with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
                r_script_path = tmp_r.name

                tmp_r.write(f"""
library(RCy3)
library(readr)

# ----------------- 연결 확인 -----------------
ping <- tryCatch(cytoscapePing(), error=function(e) NULL)
if (is.null(ping)) stop("Cytoscape not reachable at localhost:1234")

# ----------------- 콤보 리스트 불러오기 -----------------
combo_names <- readRDS(file.path("{input_root}", "{combo_file}"))

result_root <- "{input_root}"
out_dir     <- "{output_dir}"

for (nm in combo_names) {{
    combo_dir <- file.path(result_root, nm)
    f <- file.path(combo_dir, "filtered_gene_list.csv")
    if (!file.exists(f)) next

    df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    sym_col <- grep("^(Geneid|Gene_Symbol|SYMBOL)$", names(df), ignore.case=TRUE, value=TRUE)[1]
    if (is.na(sym_col)) next

    genes <- unique(na.omit(trimws(as.character(df[[sym_col]]))))
    if (length(genes) < 2) next

    gene_str <- paste(genes, collapse=",")

    # STRING network 생성
    cmd <- sprintf('string protein query query="%s" taxonID=%d cutoff=%s limit=%d',
                   gene_str, {taxon_id}, {cutoff}, {limit})
    commandsRun(cmd)

    # 네트워크 이름 변경
    net_suid <- getNetworkSuid()
    net_name <- paste0("STRING_", nm)
    renameNetwork(net_name, network=net_suid)

    # SVG 저장
    combo_out <- file.path(out_dir, nm)
    if (!dir.exists(combo_out)) dir.create(combo_out, recursive=TRUE)
    out_file <- file.path(combo_out, paste0("STRING_", nm, ".svg"))
    fitContent()
    exportImage(out_file, type="SVG")
}}
""")

            # Rscript 실행
            result = subprocess.run(
                ["Rscript", r_script_path],
                capture_output=True,
                text=True,
                encoding="utf-8"
            )
            if result.returncode == 0:
                st.success("STRING network generation completed!")
                if result.stdout:
                    st.text(result.stdout)
            else:
                st.error("R script execution failed.")
                if result.stderr:
                    st.text(result.stderr)

    # ----------------- Result -----------------
    with result_tab:
        if os.path.exists(output_dir):
            combo_dirs = [os.path.join(output_dir, d) for d in os.listdir(output_dir) 
                          if os.path.isdir(os.path.join(output_dir, d))]
            if combo_dirs:
                st.markdown("### STRING Network Results")
                for d in combo_dirs:
                    svgs = [f for f in os.listdir(d) if f.endswith(".svg")]
                    for f in svgs:
                        st.write(f"**{os.path.basename(d)}: {f}**")
                        st.image(os.path.join(d, f), use_container_width=True)
            else:
                st.info("No combo directories found.")
        else:
            st.info("Output directory does not exist.")

    # ----------------- Download -----------------
    with download_tab:
        if os.path.exists(output_dir) and os.listdir(output_dir):
            with tempfile.TemporaryDirectory() as tmpdir:
                zip_path = shutil.make_archive(os.path.join(tmpdir, "STRING_results"), "zip", output_dir)
                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download STRING results (ZIP)",
                        data=f,
                        file_name="STRING_results.zip",
                        mime="application/zip"
                    )
        else:
            st.info("No files to download.")
