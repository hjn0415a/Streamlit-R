import os
import subprocess
import tempfile
import shutil
import streamlit as st
import pandas as pd
import pyreadr

from src.common.common import page_setup
params = page_setup()

st.title("GSEA Plot")

# ----------------- Main Tabs -----------------
main_tabs = st.tabs(["📊 GSEA Plot (total)", "📊 GSEA Term Plot"])
total_tab, term_tab = main_tabs

# ----------------- 1️⃣ Total GSEA Plot (topN terms) -----------------
with total_tab:

    # ----------------- Sub Tabs -----------------
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Configure -----------------
    with configure_tab:
        input_dir = "/data/GSEA"
        output_dir = "/data/GSEA/gseaplot_total"
        os.makedirs(output_dir, exist_ok=True)

        topN = st.number_input("Top N terms to plot", value=10, step=1, min_value=1)
        width = st.number_input("Plot width", value=12.0, step=0.5)
        height = st.number_input("Plot height", value=8.0, step=0.5)

    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run Total gseaplot2 Generation"):
            with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
                r_script_path = tmp_r.name
                tmp_r.write(f"""
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

input_dir  <- "{input_dir}"
output_dir <- "{output_dir}"

files <- c(BP = "gse_BP.rds",
           CC = "gse_CC.rds",
           MF = "gse_MF.rds")

topN <- {topN}

for (ont in names(files)) {{
  rds_path <- file.path(input_dir, files[ont])
  if (!file.exists(rds_path)) next

  gse <- try(readRDS(rds_path), silent = TRUE)
  if (inherits(gse, "try-error") || !inherits(gse, "gseaResult")) next
  if (is.null(gse@result) || nrow(gse@result) == 0) next

  res <- as.data.frame(gse@result)
  res[] <- lapply(res, function(x) if (inherits(x, "Rle")) as.vector(x) else x)

  if (!all(c("ID","p.adjust") %in% names(res))) next
  res$p.adjust <- suppressWarnings(as.numeric(res$p.adjust))
  ord <- order(res$p.adjust, na.last = NA)
  if (!length(ord)) next
  sel <- ord[ seq_len(min(topN, length(ord))) ]
  ids <- res$ID[sel]

  p <- gseaplot2(
    gse,
    geneSetID    = ids,
    pvalue_table = TRUE,
    title        = sprintf("Top %d enriched GO:%s terms", length(ids), ont)
  )

  ggsave(
    filename = file.path(output_dir, sprintf("gseaplot2_%s_top%d.svg", ont, length(ids))),
    plot     = p,
    width    = {width},
    height   = {height}
  )
}}
""")
            result = subprocess.run(
                ["Rscript", r_script_path],
                capture_output=True,
                text=True,
                encoding="utf-8"
            )
            if result.returncode == 0:
                st.success("Total gseaplot2 generation completed!")
                if result.stdout:
                    st.text(result.stdout)
            else:
                st.error("R script execution failed.")
                if result.stderr:
                    st.text(result.stderr)

    # ----------------- Result -----------------
    with result_tab:
        if os.path.exists(output_dir):
            ontology_tabs = st.tabs(["BP", "CC", "MF"])
            for ont_tab, ont in zip(ontology_tabs, ["BP", "CC", "MF"]):
                with ont_tab:
                    svgs = [f for f in os.listdir(output_dir) if f.endswith(".svg") and f"_{ont}_" in f]
                    if svgs:
                        for f in svgs:
                            st.write(f"**{f}**")
                            st.image(os.path.join(output_dir, f), width=850)
                    else:
                        st.info(f"No SVG plots found for {ont}.")
        else:
            st.info("Output directory does not exist.")

    # ----------------- Download -----------------
    with download_tab:
        if os.path.exists(output_dir) and os.listdir(output_dir):
            with tempfile.TemporaryDirectory() as tmpdir:
                zip_path = shutil.make_archive(os.path.join(tmpdir, "gseaplot2_total_results"), "zip", output_dir)
                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download Total gseaplot2 Results (ZIP)",
                        data=f,
                        file_name="gseaplot2_total_results.zip",
                        mime="application/zip"
                    )
        else:
            st.info("No files to download.")

# ----------------- 2️⃣ GSEA Term Plot (per selected term) -----------------
with term_tab:
    tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])

    # ----------------- Configure -----------------
    with tabs[0]:
        st.subheader("Configure Parameters")
        input_dir = "/data/GSEA"
        output_dir = "/data/GSEA/gseaplot"
        os.makedirs(output_dir, exist_ok=True)

        ont = st.selectbox("Select ontology", ["BP", "CC", "MF"], index=0)
        # CSV 파일 경로 설정
        ont_files = {"BP": "gse_BP.csv", "CC": "gse_CC.csv", "MF": "gse_MF.csv"}
        csv_path = os.path.join(input_dir, ont_files[ont])
    
        if os.path.exists(csv_path):
            try:
                df = pd.read_csv(csv_path)
                if df.empty:
                    st.info(f"{ont} CSV 파일이 비어 있습니다.")
                else:
                    st.write(f"Showing results for **{ont}**")
                    # 1-based 인덱스로 표시
                    df_display = df.copy()
                    df_display.index = range(1, len(df_display) + 1)
                    st.dataframe(df_display)

                    idx = st.number_input(
                            "Row index (1-based) for GSEA Term Plot",
                            min_value = 1,
                            max_value = len(df_display),
                            value = 1,
                            step = 1
                            )
                    st.session_state["selected_idx"] = idx
            except Exception as e:
                st.error(f"CSV 파일 읽기 실패: {e}")
        else:
            st.warning(f"파일이 존재하지 않습니다: {csv_path}")

    # ----------------- Run -----------------
    with tabs[1]:
        st.subheader("Run GSEA Term Plot")
        idx = st.session_state.get("selected_idx", 1)

        if st.button("Run plot_gsea_term"):
            with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
                r_script_path = tmp_r.name
                tmp_r.write(f"""
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(cowplot)

input_dir  <- "{input_dir}"
output_dir <- "{output_dir}"

plot_gsea_term <- function(ont = "BP", idx = 1,
                           input_dir, output_dir,
                           width = 8, height = 8) {{
  ont_files <- c(BP = "gse_BP.rds", CC = "gse_CC.rds", MF = "gse_MF.rds")
  if (!ont %in% names(ont_files)) stop("ont는 'BP', 'CC', 'MF' 중 하나여야 합니다.")
  
  rds_path <- file.path(input_dir, ont_files[[ont]])
  if (!file.exists(rds_path)) stop("파일이 존재하지 않습니다: ", rds_path)
  
  gse <- readRDS(rds_path)
  if (!inherits(gse, "gseaResult")) stop("읽은 객체가 gseaResult가 아닙니다.")
  
  res <- as.data.frame(gse)
  if (nrow(res) == 0) stop("결과가 비어 있습니다.")
  if (idx < 1 || idx > nrow(res)) stop(sprintf("idx 범위는 1 ~ %d 입니다.", nrow(res)))
  
  term_id   <- res$ID[idx]
  term_desc <- res$Description[idx]
  
  p <- gseaplot2(gse, geneSetID = term_id, title = term_desc)
  
  out_name <- sprintf("gseaplot_%s_idx%d_%s.svg", ont, idx, gsub("[:/\\\\]+", "_", term_id))
  ggsave(file.path(output_dir, out_name), plot = p, width = width, height = height)
  
  invisible(p)
}}

plot_gsea_term("{ont}", {idx}, input_dir, output_dir)
""")
            result = subprocess.run(["Rscript", r_script_path], capture_output=True, text=True, encoding="utf-8")
            if result.returncode == 0:
                st.success("plot_gsea_term 실행 완료!")
            else:
                st.error("실행 중 오류 발생")
                st.text(result.stderr)

    # ----------------- Result -----------------
    with tabs[2]:
        st.subheader("Generated Plots")
        if os.path.exists(output_dir):
            images = [f for f in os.listdir(output_dir) if f.endswith(".svg")]
            if images:
                for img_file in images:
                    st.markdown(f"**{img_file}**")
                    st.image(os.path.join(output_dir, img_file), width=1000)
            else:
                st.info("아직 생성된 플롯이 없습니다.")
        else:
            st.warning("Output directory does not exist.")

    # ----------------- Download -----------------
    with tabs[3]:
        st.subheader("Download Plots")
        if os.path.exists(output_dir):
            files = [f for f in os.listdir(output_dir) if f.endswith(".svg")]
            if files:
                for f in files:
                    file_path = os.path.join(output_dir, f)
                    with open(file_path, "rb") as fp:
                        st.download_button(
                            label=f"Download {f}",
                            data=fp,
                            file_name=f,
                            mime="image/svg+xml"
                        )
            else:
                st.info("다운로드할 플롯이 없습니다.")
        else:
            st.warning("Output directory does not exist.")
