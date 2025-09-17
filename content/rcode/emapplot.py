import os
import subprocess
import tempfile
import streamlit as st

from src.common.common import page_setup

params = page_setup()

st.markdown("## 📊 Enrichment Map (emap)")

base_dir = "/data"
figure_dir = os.path.join(base_dir, "figure")
os.makedirs(figure_dir, exist_ok=True)

with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
    r_script_path = tmp_r.name
    tmp_r.write(f"""
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(readr)
    library(ggplot2)
    library(svglite)

    a <- 1

    base_dir   <- "{base_dir}"
    figure_dir <- file.path(base_dir, "figure")
    dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

    # 1) 데이터 불러오기
    gene_data <- read_csv(file.path(base_dir, "sample_data.csv"))

    # 2) Gene name, FC2 추출
    gene_symbols <- unique(trimws(unlist(strsplit(as.character(gene_data$`Gene names`), "\\\\s*;\\\\s*"))))
    gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]

    # 3) SYMBOL → ENTREZID 변환
    converted <- bitr(gene_symbols,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Hs.eg.db)
    entrez_ids <- unique(na.omit(converted$ENTREZID))

    # 5) emapplot 함수
    make_emap <- function(genes, ont, out_svg, show_n = 20) {{
      ego <- enrichGO(gene          = genes,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = ont,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE)
      if (is.null(ego) || nrow(ego@result) < 2) return(NULL)
      ego <- pairwise_termsim(ego)
      p <- emapplot(ego, showCategory = min(show_n, nrow(ego@result)), layout = "kk")
      ggsave(out_svg, p, width = 9, height = 7, device = svglite::svglite)
    }}

    make_emap(entrez_ids, "BP", file.path(figure_dir,"emap_BP.svg"))
    make_emap(entrez_ids, "CC", file.path(figure_dir,"emap_CC.svg"))
    make_emap(entrez_ids, "MF", file.path(figure_dir,"emap_MF.svg"))
    """)

try:
    subprocess.run(["Rscript", r_script_path], check=True)
    st.success("enrichment plots generated!")

    # 결과 SVG 보여주기
    for fname in ["emap_BP.svg", "emap_CC.svg", "emap_MF.svg"]:
        svg_file = os.path.join(figure_dir, fname)
        if os.path.exists(svg_file):
            st.image(svg_file, use_column_width=True)
except subprocess.CalledProcessError as e:
    st.error(f"R script failed: {e}")
finally:
    os.remove(r_script_path)
