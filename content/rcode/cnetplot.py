import os
import subprocess
import tempfile
import shutil
import streamlit as st
import pandas as pd

from src.common.common import page_setup
params = page_setup()

st.title("Cnetplot")

# ----------------- Main Tabs -----------------
main_tabs = st.tabs(["🧬 Cnetplot"])
with main_tabs[0]:

    # ----------------- Sub Tabs -----------------
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Configure -----------------
    with configure_tab:
        result_root = "/data/Enrichment/out"
        figure_root = "/data/Cnet"
        combo_root = "/data/Deg"
        os.makedirs(figure_root, exist_ok=True)

        # ----------------- 1️⃣ combo_names.csv 읽기 & FC/P-value 선택 -----------------
        combo_csv = os.path.join(combo_root, "combo_names.csv")
        if os.path.exists(combo_csv):
            combo_df = pd.read_csv(combo_csv)
            fc_values = sorted(list({float(c.split("_")[0][2:]) for c in combo_df["combo"]}))
            pval_values = sorted(list({float(c.split("_")[1][1:]) for c in combo_df["combo"]}))

            fc_threshold = st.selectbox("Select FC threshold", options=fc_values, index=0)
            pval_threshold = st.selectbox("Select P-value threshold", options=pval_values, index=0)
        else:
            fc_threshold = None
            pval_threshold = None

        # ----------------- 2️⃣ 나머지 파라미터 입력 -----------------
        showCategory = st.number_input(
            "Number of categories to show (show_n)", value=10, step=1
        )
        plot_width = st.number_input("Plot width", value=9.0, step=0.5)
        plot_height = st.number_input("Plot height", value=7.0, step=0.5)
        
    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run Cnetplot Generation"):
            # combo_names.csv 읽어서 선택된 FC/P-value 조합 필터링
            combo_csv = os.path.join(combo_root, "combo_names.csv")
            if os.path.exists(combo_csv):
                combo_df = pd.read_csv(combo_csv)
                selected_combos = [
                    c
                    for c in combo_df["combo"]
                    if float(c.split("_")[0][2:]) == fc_threshold
                    and float(c.split("_")[1][1:]) == pval_threshold
                ]
            else:
                selected_combos = []

            if not selected_combos:
                st.warning("No combos match the selected FC/P-value thresholds.")
            else:
                with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
                    r_script_path = tmp_r.name
                    combos_r = 'c(' + ','.join([f'"{c}"' for c in selected_combos]) + ')'

                    tmp_r.write(f"""
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(svglite)

.find_ego_rds <- function(combo_dir, ont) {{
  direct <- file.path(combo_dir, sprintf("GO_%s_ego.rds", ont))
  if (file.exists(direct)) return(direct)
  cand <- list.files(combo_dir,
                     pattern = paste0("^GO_", ont, "_ego\\\\.rds$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(cand) >= 1) return(cand[1])
  NA_character_
}}

make_cnet_from_rds_by_combo <- function(result_root,
                                        figure_root,
                                        combo_names,
                                        onts     = c("BP","CC","MF"),
                                        show_n   = {showCategory},
                                        width    = {plot_width},
                                        height   = {plot_height},
                                        fc_vec   = NULL,
                                        circular = TRUE,
                                        layout   = "kk") {{
  if (!dir.exists(figure_root)) dir.create(figure_root, recursive = TRUE)

  for (nm in combo_names) {{
    combo_dir <- file.path(result_root, nm)
    out_dir   <- file.path(figure_root, nm)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    for (ont in onts) {{
      rds_path <- .find_ego_rds(combo_dir, ont)
      if (is.na(rds_path)) next

      ego <- readRDS(rds_path)
      if (is.null(ego) || is.null(ego@result) || nrow(ego@result) < 1) next

      k <- min(show_n, nrow(ego@result))
      p <- if (is.null(fc_vec)) {{
        cnetplot(ego, showCategory = k, circular = circular, layout = layout)
      }} else {{
        cnetplot(ego, showCategory = k, circular = circular, layout = layout,
                 foldChange = fc_vec)
      }}

      ggsave(file.path(out_dir, sprintf("cnet_%s.svg", ont)),
             p, width = width, height = height, device = svglite::svglite)
    }}
  }}
}}

combo_names <- {combos_r}

make_cnet_from_rds_by_combo(
  result_root = "{result_root}",
  figure_root = "{figure_root}",
  combo_names = combo_names
)
""")
                result = subprocess.run(
                    ["Rscript", r_script_path],
                    capture_output=True,
                    text=True,
                    encoding="utf-8"
                )
                if result.returncode == 0:
                    st.success("Cnetplot generation completed!")
                    st.text(result.stdout)
                else:
                    st.error("R script execution failed.")
                    st.text(result.stderr)

        # ----------------- Result -----------------
        with result_tab:
            combo_csv = os.path.join(combo_root, "combo_names.csv")
            if os.path.exists(combo_csv):
                combos = pd.read_csv(combo_csv)["combo"].tolist()
                selected_combos = [
                    c
                    for c in combos
                    if float(c.split("_")[0][2:]) == fc_threshold
                    and float(c.split("_")[1][1:]) == pval_threshold
                ]

                if selected_combos:
                    ontology_tabs = st.tabs(["BP", "CC", "MF"])
                    for ont_tab, ont in zip(ontology_tabs, ["BP", "CC", "MF"]):
                        with ont_tab:
                            st.subheader(f"Ontology: {ont}")
    
                            # combo별 2개씩 가로 배치
                            for i in range(0, len(selected_combos), 2):
                                combo_pair = selected_combos[i:i+2]
                                cols = st.columns(len(combo_pair))

                                for col, combo in zip(cols, combo_pair):
                                    with col:
                                        st.markdown(f"**{combo}**")
                                        plot_file = os.path.join(figure_root, combo, f"cnet_{ont}.svg")
                                        if os.path.exists(plot_file):
                                            st.image(plot_file, use_container_width=True)
                                        else:
                                            st.info(f"No {ont} plot found for {combo}")
                else:
                    st.info("No combos found for selected thresholds.")
            else:
                st.warning("⚠️ combo_names.csv not found.")

    # ----------------- Download -----------------
    with download_tab:
        combo_csv = os.path.join(combo_root, "combo_names.csv")
        if os.path.exists(combo_csv):
            combos = pd.read_csv(combo_csv)["combo"].tolist()
            selected_combos = [
                c
                for c in combos
                if float(c.split("_")[0][2:]) == fc_threshold
                and float(c.split("_")[1][1:]) == pval_threshold
            ]
            if selected_combos:
                with tempfile.TemporaryDirectory() as tmpdir:
                    for combo in selected_combos:
                        src = os.path.join(figure_root, combo)
                        dst = os.path.join(tmpdir, combo)
                        if os.path.exists(src):
                            shutil.copytree(src, dst)
                    zip_path = shutil.make_archive(os.path.join(tmpdir, "Cnetplots_combos"), "zip", tmpdir)
                    with open(zip_path, "rb") as f:
                        st.download_button(
                            label="Download Cnetplot Results (ZIP)",
                            data=f,
                            file_name="Cnetplots_combos.zip",
                            mime="application/zip"
                        )
