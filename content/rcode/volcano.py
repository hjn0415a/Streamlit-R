import os
import subprocess
import tempfile
import streamlit as st

from src.common.common import page_setup

params = page_setup()

st.markdown("## Volcano Plot Example using R")

fc_cutoff_volcano = st.number_input("Fold Change Cutoff (log2)", value=1.0, key="volcano_fc_cutoff")
pval_cutoff_volcano = st.number_input("P-value Cutoff", value=0.05, format="%.4f", key="volcano_pval_cutoff")

excel_path = "/data/example_data.csv"  # CSV 파일 경로
output_svg = excel_path.replace(".csv", "_volcano.svg")

with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
    r_script_path = tmp_r.name
    tmp_r.write(f"""
library(readr)
library(ggplot2)

data <- read_csv('{excel_path}')

# NA 제거 및 log2FC 계산
data <- data[!is.na(data$foldchange) & data$foldchange > 0, ]
data$log2FC <- log2(data$foldchange)

fc_cutoff <- {fc_cutoff_volcano}
pval_cutoff <- {pval_cutoff_volcano}
p_cutoff <- -log10(pval_cutoff)

# 그룹 지정
data$group <- "NS"
data$group[data$log2FC > fc_cutoff & -log10(data$pvalue) > p_cutoff] <- "Up"
data$group[data$log2FC < -fc_cutoff & -log10(data$pvalue) > p_cutoff] <- "Down"

color_palette <- c("NS" = "gray", "Up" = "red", "Down" = "blue")

volcano_plot <- ggplot(data, aes(x = log2FC, y = -log10(pvalue), color = group)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = color_palette) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = p_cutoff, linetype = "dashed", color = "gray") +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 pvalue") +
    theme_minimal()

ggsave(filename='{output_svg}', plot=volcano_plot, width=8, height=6, dpi=300, device='svg')
    """)

subprocess.run(["Rscript", r_script_path], check=True)
st.success("Volcano plot generated!")
st.image(output_svg, use_column_width=True)
os.remove(r_script_path)

# ----------------- EnhancedVolcano Plot -----------------
st.markdown("## EnhancedVolcano Plot Example using R")

fc_cutoff_enhanced = st.number_input("Fold Change Cutoff (log2)", value=1.0, key="enhanced_fc_cutoff")
pval_cutoff_enhanced = st.number_input("P-value Cutoff", value=0.05, format="%.4f", key="enhanced_pval_cutoff")

output_svg_enhanced = excel_path.replace(".csv", "_enhanced_volcano.svg")

with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tmp_r:
    r_script_path = tmp_r.name
    tmp_r.write(f"""
library(readr)
library(EnhancedVolcano)

data <- read_csv('{excel_path}')

# NA 제거 및 log2FC 계산
data <- data[!is.na(data$foldchange) & data$foldchange > 0, ]
data$log2FC <- log2(data$foldchange)

# EnhancedVolcano용 데이터 구성
res <- data.frame(
    log2FoldChange = data$log2FC,
    pvalue = data$pvalue
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
    title = 'Enhanced Plot',
    subtitle = 'Sample',
    pointSize = 2.5,
    labSize = 2.5
)

dev.off()
    """)

subprocess.run(["Rscript", r_script_path], check=True)
st.success("EnhancedVolcano plot generated!")
st.image(output_svg_enhanced, use_column_width=True)
os.remove(r_script_path)
