import os
import subprocess
import tempfile
import shutil
import streamlit as st

from src.common.common import page_setup
params = page_setup()

st.title("Heatmaplike functional classification")

# ----------------- Main Tabs -----------------
main_tabs = st.tabs(["📊 GSEA Heatplot"])
with main_tabs[0]:

    # ----------------- Sub Tabs -----------------
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Configure -----------------
    with configure_tab:
        # 경로 고정
        edox_dir = "/data/GSEA"          # gse_BP.rds, gse_CC.rds, gse_MF.rds
        csv_path = "/data/example_data.csv"
        output_dir = "/data/GSEA/heatplot"
        os.makedirs(output_dir, exist_ok=True)

        # 사용자 입력 파라미터
        top_pathways = st.number_input("Max pathways to display", value=5, step=1, min_value=1)
        top_genes_per_pathway = st.number_input("Max genes per pathway", value=20, step=1, min_value=1)
        width = st.number_input("Plot width", value=12.0, step=0.5)
        height = st.number_input("Plot height", value=6.0, step=0.5)

        # 고정값
        max_setsize = 50

    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run Heatplot Generation"):
            with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
                r_script_path = tmp_r.name
                tmp_r.write(f"""
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(cowplot)
library(org.Hs.eg.db)

# 사용자 파라미터
csv_path <- "{csv_path}"
edox_dir <- "{edox_dir}"
out_dir  <- "{output_dir}"
top_pathways <- {top_pathways}
top_genes <- {top_genes_per_pathway}
max_setsize <- {max_setsize}

# 1) CSV에서 foldchange 불러오기
df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("Geneid" %in% names(df), "foldchange" %in% names(df))

fc_vec <- setNames(log2(df$foldchange + 1e-8), df$Geneid)
fc_vec <- fc_vec[is.finite(fc_vec)]

# ID 자동 맞춤 함수
harmonize_fc <- function(edox, fc_named) {{
  core_ids <- unique(unlist(strsplit(edox@result$core_enrichment, "/")))
  if (length(core_ids) == 0 || all(is.na(core_ids))) return(fc_named)
  if (sum(names(fc_named) %in% core_ids) > 0) return(fc_named)
  
  is_entrez <- all(grepl("^[0-9]+$", head(core_ids[!is.na(core_ids)], 50)))
  if (is_entrez) {{
    mp <- suppressMessages(bitr(names(fc_named), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    if (!nrow(mp)) return(numeric(0))
    mp <- mp[!duplicated(mp$SYMBOL), ]
    out <- fc_named[mp$SYMBOL]; names(out) <- mp$ENTREZID
    out[!is.na(names(out))]
  }} else {{
    mp <- suppressMessages(bitr(names(fc_named), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db))
    if (!nrow(mp)) return(numeric(0))
    mp <- mp[!duplicated(mp$ENTREZID), ]
    out <- fc_named[mp$ENTREZID]; names(out) <- mp$SYMBOL
    out[!is.na(names(out))]
  }}
}}

# 2) gseGO heatplot
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
files <- c(BP="gse_BP.rds", CC="gse_CC.rds", MF="gse_MF.rds")

for (ont in names(files)) {{
  rds <- file.path(edox_dir, files[ont])
  if (!file.exists(rds)) next
  edox <- try(readRDS(rds), silent = TRUE)
  if (inherits(edox, "try-error") || nrow(edox@result) == 0) next
  
  if (!is.null(max_setsize) && "setSize" %in% names(edox@result)) {{
    edox@result <- edox@result[edox@result$setSize <= max_setsize, , drop = FALSE]
    if (nrow(edox@result) == 0) next
  }}
  
  fc_use <- harmonize_fc(edox, fc_vec)
  if (!length(fc_use)) next
  
  fc_use <- fc_use[is.finite(fc_use)]
  res_top <- head(edox@result[order(edox@result$p.adjust), , drop = FALSE], top_pathways)
  core_top <- unique(unlist(strsplit(res_top$core_enrichment, "/")))
  fc_use <- fc_use[names(fc_use) %in% core_top]
  fc_use <- fc_use[order(abs(fc_use), decreasing = TRUE)]
  fc_use <- head(fc_use, top_genes)
  
  p1 <- heatplot(edox, showCategory = top_pathways) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  p2 <- heatplot(edox, foldChange = fc_use, showCategory = top_pathways) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  g <- cowplot::plot_grid(p1, p2, ncol = 1, labels = c("A","B"))
  ggsave(file.path(out_dir, sprintf("heatplot_%s_top%dgenes.svg", ont, top_genes)), g,
         width = {width}, height = {height})
}}
""")

            result = subprocess.run(
                ["Rscript", r_script_path],
                capture_output=True,
                text=True,
                encoding="utf-8"
            )
            if result.returncode == 0:
                st.success("Heatplot generation completed!")
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
                    st.subheader(f"Ontology: {ont}")
                    ont_svgs = [f for f in os.listdir(output_dir) if f.endswith(".svg") and ont in f]
                    if ont_svgs:
                        for i in range(0, len(ont_svgs), 2):
                            svg_pair = ont_svgs[i:i+2]
                            cols = st.columns(len(svg_pair))
                            for col, svg_file in zip(cols, svg_pair):
                                with col:
                                    st.markdown(f"**{svg_file}**")
                                    st.image(os.path.join(output_dir, svg_file), use_container_width=True)
                    else:
                        st.info(f"No {ont} heatplots found.")
        else:
            st.info("Output directory does not exist.")

    # ----------------- Download -----------------
    with download_tab:
        if os.path.exists(output_dir) and os.listdir(output_dir):
            with tempfile.TemporaryDirectory() as tmpdir:
                zip_path = shutil.make_archive(os.path.join(tmpdir, "heatplot_results"), "zip", output_dir)
                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download Heatplot Results (ZIP)",
                        data=f,
                        file_name="heatplot_results.zip",
                        mime="application/zip"
                    )
        else:
            st.info("No files to download.")
