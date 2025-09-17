import os
import subprocess
import tempfile
import streamlit as st

from src.common.common import page_setup

params = page_setup()

st.markdown("## 🕸️ Cnet Plot (cnet)")

excel_path_cnet = "/data/data3_fc2_raw.p.xlsx"
output_dir_cnet = "/data/go_results"
figure_dir_cnet = os.path.join(output_dir_cnet, "cnet_figure")
os.makedirs(figure_dir_cnet, exist_ok=True)

with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
    r_script_path = tmp_r.name
    tmp_r.write(f"""
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(enrichplot)
    library(readxl)
    library(ggplot2)
    library(svglite)

    excel_path <- "{excel_path_cnet}"
    result_dir <- "{output_dir_cnet}"
    figure_dir <- "{figure_dir_cnet}"
    dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

    # 데이터 읽기
    gene_data <- read_excel(excel_path, sheet = "data3_fc2_raw.p")

    # log2FoldChange 계산
    gene_data$log2FC <- log2((gene_data$Dd_FPKM + 1) / (gene_data$Db_FPKM + 1))

    symbols_raw <- gene_data$Gene_Symbol
    symbols_vec <- unique(trimws(unlist(strsplit(as.character(symbols_raw), "\\\\s*;\\\\s*"))))
    symbols_vec <- symbols_vec[symbols_vec != ""]

    fc_vec <- gene_data$log2FC
    names(fc_vec) <- gene_data$Gene_Symbol

    converted <- bitr(symbols_vec, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
    entrez_ids <- unique(na.omit(converted$ENTREZID))

    # cnetplot 함수 정의
    make_cnet_svg <- function(genes, ont, out_csv, out_svg, show_n=10, fc_vec=NULL) {{
        ego <- enrichGO(gene=genes, OrgDb=org.Mm.eg.db, keyType="ENTREZID",
                        ont=ont, pAdjustMethod="BH", pvalueCutoff=0.9, readable=TRUE)
        if (is.null(ego) || nrow(ego@result) < 2) {{
            write.csv(data.frame(), out_csv, row.names=FALSE)
            return(invisible(NULL))
        }}
        write.csv(ego@result, out_csv, row.names=FALSE)
        p <- cnetplot(ego, showCategory=min(show_n, nrow(ego@result)), foldChange=fc_vec,
                      circular=TRUE, layout="kk") +
             scale_color_gradient2(name="log2FC", low="steelblue", mid="white", high="firebrick", midpoint=0)
        ggsave(out_svg, p, width=9, height=7, device=svglite::svglite)
    }}

    # BP / CC / MF
    make_cnet_svg(entrez_ids, "BP",
                  file.path(result_dir, "GO_BP_result_2.csv"),
                  file.path(figure_dir, "cnetplot_BP.svg"),
                  show_n=10, fc_vec=fc_vec)
    make_cnet_svg(entrez_ids, "CC",
                  file.path(result_dir, "GO_CC_result_2.csv"),
                  file.path(figure_dir, "cnetplot_CC.svg"),
                  show_n=10, fc_vec=fc_vec)
    make_cnet_svg(entrez_ids, "MF",
                  file.path(result_dir, "GO_MF_result_2.csv"),
                  file.path(figure_dir, "cnetplot_MF.svg"),
                  show_n=10, fc_vec=fc_vec)
    """)

subprocess.run(["Rscript", r_script_path], check=True)
st.success("cnet plots generated!")
st.image(os.path.join(figure_dir_cnet, "cnetplot_BP.svg"), use_column_width=True)
