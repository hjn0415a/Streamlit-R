import streamlit as st
import subprocess
import tempfile
import os

from src.common.common import page_setup

params = page_setup()

tabs = st.tabs(["Volcano Plots", "Enrichment", "Heatmap", "PCA", "Emap", "Cnet"])

with tabs[0]:
    # ----------------- Volcano Plot -----------------
    st.markdown("## Volcano Plot Example using R")

    fc_cutoff_volcano = st.number_input("Fold Change Cutoff (log2)", value=1.0, key="volcano_fc_cutoff")
    pval_cutoff_volcano = st.number_input("P-value Cutoff", value=0.05, format="%.4f", key="volcano_pval_cutoff")

    excel_path = "/data/data2_All_data.xlsx"
    sheet_name = "data2"

    output_svg = excel_path.replace(".xlsx", "_volcano.svg")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
        r_script_path = tmp_r.name
        tmp_r.write(f"""
    library(readxl)
    library(ggplot2)

    data <- read_excel('{excel_path}', sheet = '{sheet_name}')

    data <- data[!is.na(data$Dd_FPKM) & !is.na(data$Db_FPKM) & data$Db_FPKM > 0, ]
    data$log2FC <- log2(data$Dd_FPKM / data$Db_FPKM)

    fc_cutoff <- {fc_cutoff_volcano}
    pval_cutoff <- {pval_cutoff_volcano}
    p_cutoff <- -log10({pval_cutoff_volcano})

    data$group <- "NS"
    data$group[data$log2FC > fc_cutoff & -log10(data$`Dd/Db.raw.pval`) > p_cutoff] <- "Up"
    data$group[data$log2FC < -fc_cutoff & -log10(data$`Dd/Db.raw.pval`) > p_cutoff] <- "Down"

    color_palette <- c("NS" = "gray", "Up" = "red", "Down" = "blue")

    volcano_plot <- ggplot(data, aes(x = log2FC, y = -log10(`Dd/Db.raw.pval`), color = group)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = color_palette) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = p_cutoff, linetype = "dashed", color = "gray") +
    labs(title = "Volcano Plot (FDR based)", x = "log2 Fold Change (N_Dd / N_Db)", y = "-log10 pval") +
    theme_minimal()

    ggsave(filename='{output_svg}', plot=volcano_plot, width=4, height=3, dpi=300, device='svg')
                    """)
        
    subprocess.run(["Rscript", r_script_path], check=True)

    st.success("Volcano plot generated!")
    st.image(output_svg, use_column_width=True)

    os.remove(r_script_path)

    # ----------------- EnhancedVolcano Plot -----------------
    st.markdown("## EnhancedVolcano Plot Example using R")

    fc_cutoff_enhanced = st.number_input(
        "Fold Change Cutoff (log2)", value=1.0, key="enhanced_fc_cutoff"
    )
    pval_cutoff_enhanced = st.number_input(
        "P-value Cutoff", value=0.05, format="%.4f", key="enhanced_pval_cutoff"
    )

    output_svg_enhanced = excel_path.replace(".xlsx", "_enhanced_volcano.svg")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
        r_script_path = tmp_r.name
        tmp_r.write(f"""
    library(readxl)
    library(EnhancedVolcano)

    data <- read_excel('{excel_path}', sheet = '{sheet_name}')
    data <- data[!is.na(data$Dd_FPKM) & !is.na(data$Db_FPKM) & data$Db_FPKM > 0, ]
    data$log2FoldChange <- log2(data$Dd_FPKM / data$Db_FPKM)

    res <- data.frame(
    log2FoldChange = data$log2FoldChange,
    pvalue = data$`Dd/Db.raw.pval`
    )
    rownames(res) <- data$Gene_Symbol
    res <- res[!is.na(res$log2FoldChange) & is.finite(res$log2FoldChange), ]

    svg('{output_svg_enhanced}', width=10, height=8)

    EnhancedVolcano(res,
        lab = NA,
        x = 'log2FoldChange',
        y = 'pvalue',
        pCutoff = {pval_cutoff_enhanced},
        FCcutoff = {fc_cutoff_enhanced},
        title = 'Volcano Plot (raw p-value)',
        subtitle = 'Comparison using N_Dd vs N_Db',
        pointSize = 2.5,
        labSize = 2.5
    )

    dev.off()
    """)

    subprocess.run(["Rscript", r_script_path], check=True)
    st.success("EnhancedVolcano plot generated!")
    st.image(output_svg_enhanced, use_column_width=True)
    os.remove(r_script_path)

with tabs[1]:
    # ----------------- Go Enrichment Analysis -----------------
    st.markdown("## GO Enrichment Analysis (ClusterProfiler)")

    excel_path_go = "/data/data3_fc2_raw.p.xlsx"
    output_dir_go = "/data/go_results"
    figure_dir_go = os.path.join(output_dir_go, "figure")
    os.makedirs(figure_dir_go, exist_ok=True)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
        r_script_path = tmp_r.name
        tmp_r.write(f"""
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(readxl)
library(ggplot2)

gene_data <- read_excel('{excel_path_go}', sheet = 'data3_fc2_raw.p')
gene_symbols <- unique(gene_data$Gene_Symbol)

converted <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
write.csv(converted, file = file.path('{output_dir_go}', "converted_ids.csv"), row.names = FALSE)

entrez_ids <- converted$ENTREZID

ego_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.9, readable = TRUE)
ego_cc <- enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.9, readable = TRUE)
ego_mf <- enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.9, readable = TRUE)

svg(file.path('{figure_dir_go}', "GO_BP.svg"), width = 8, height = 6)
dotplot(ego_bp, showCategory = 10, x = "GeneRatio", color = "p.adjust") + ggtitle("GO BP")
dev.off()

svg(file.path('{figure_dir_go}', "GO_CC.svg"), width = 8, height = 6)
dotplot(ego_cc, showCategory = 10, x = "GeneRatio", color = "p.adjust") + ggtitle("GO CC")
dev.off()

svg(file.path('{figure_dir_go}', "GO_MF.svg"), width = 8, height = 6)
dotplot(ego_mf, showCategory = 10, x = "GeneRatio", color = "p.adjust") + ggtitle("GO MF")
dev.off()

write.csv(ego_bp@result, file = file.path('{output_dir_go}', "GO_BP_result.csv"), row.names = FALSE)
write.csv(ego_cc@result, file = file.path('{output_dir_go}', "GO_CC_result.csv"), row.names = FALSE)
write.csv(ego_mf@result, file = file.path('{output_dir_go}', "GO_MF_result.csv"), row.names = FALSE)
        """)
    result = subprocess.run(["Rscript", r_script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE
                            
                            , text=True)
    st.text(result.stdout)
    st.text(result.stderr)

    st.success("GO enrichment analysis completed!")

with tabs[2]:
    # ----------------- Heatmap -----------------
    st.markdown("## Heatmap (pheatmap)")

    excel_path_heatmap = "/data/HEATMAP_significant.xlsx"  # 컨테이너 안에 둘 파일 경로
    output_svg_heatmap = "/data/heatmap2.svg"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
        r_script_path = tmp_r.name
        tmp_r.write(f"""
library(pheatmap)
library(readxl)
library(svglite)

data <- read_excel("{excel_path_heatmap}")

gene_names <- data[[1]]
data <- data[ , -1]
mat <- as.matrix(data)
rownames(mat) <- gene_names

annotation_col <- data.frame(
  Group = factor(c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"))
)
rownames(annotation_col) <- colnames(mat)

svglite("{output_svg_heatmap}", width = 8, height = 10)

pheatmap(mat,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue", "white", "red"))(100)
)

dev.off()
        """)

    result = subprocess.run(
        ["Rscript", r_script_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    st.text(result.stdout)
    st.text(result.stderr)

    if os.path.exists(output_svg_heatmap):
        st.success("Heatmap generated successfully!")
        st.image(output_svg_heatmap, use_column_width=True)
    else:
        st.error("Heatmap generation failed. Check logs above.")

    os.remove(r_script_path)

with tabs[3]:
    # ----------------- PCA Plot -----------------
    st.markdown("## PCA Plot (factoextra + ggrepel)")

    excel_path_pca = "/data/proteinGroups.xlsx"  # 컨테이너 안 경로
    output_svg_pca = "/data/PCA_plot6.svg"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
        r_script_path = tmp_r.name
        tmp_r.write(f"""
library(readxl)
library(factoextra)
library(ggrepel)
library(svglite)

dat <- as.data.frame(read_excel("{excel_path_pca}"))
lfq <- dat[, -1]   # 첫 번째 열 Protein ID 제외

X <- t(as.matrix(lfq))
rownames(X) <- colnames(lfq)

Xz <- scale(X, center = TRUE, scale = TRUE)
Xz <- Xz[, colSums(is.na(Xz)) == 0, drop = FALSE]

pca_res <- prcomp(Xz, center = FALSE, scale. = FALSE)
rownames(pca_res$x) <- rownames(X)

groups <- factor(c(rep("Group1", 3), rep("Group2", 3)))
df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                 sample = rownames(pca_res$x), group = groups)

svglite("{output_svg_pca}", width = 8, height = 6)
p <- fviz_pca_ind(
  pca_res,
  geom.ind   = "point",
  col.ind    = groups,
  pointshape = 16,
  pointsize  = 3.5,
  mean.point = FALSE,
  addEllipses= FALSE
) +
  geom_text_repel(data = df, aes(PC1, PC2, label = sample, color = group),
                  size = 4, show.legend = FALSE)
print(p)
dev.off()
        """)

    result = subprocess.run(
        ["Rscript", r_script_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    st.text(result.stdout)
    st.text(result.stderr)

    if os.path.exists(output_svg_pca):
        st.success("PCA plot generated successfully!")
        st.image(output_svg_pca, use_column_width=True)
    else:
        st.error("PCA plot generation failed. Check logs above.")

    os.remove(r_script_path)

with tabs[4]:
        # ----------------- Emap -----------------
    st.markdown("## 📊 Enrichment Map (emap)")

    base_dir = "/data"
    figure_dir = os.path.join(base_dir, "figure")
    os.makedirs(figure_dir, exist_ok=True)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
        r_script_path = tmp_r.name
        tmp_r.write(f"""
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(enrichplot)
        library(readr)
        library(ggplot2)
        library(svglite)

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

with tabs[5]:
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
