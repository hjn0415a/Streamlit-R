import os
import subprocess
import tempfile
import streamlit as st

def run_r_script(r_code: str, output_svg: str):
    """R 스크립트를 실행하고 결과 파일 생성"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
        r_script_path = tmp_r.name
        tmp_r.write(r_code)
    subprocess.run(["Rscript", r_script_path], check=True)
    os.remove(r_script_path)
    return os.path.exists(output_svg)

st.title("Volcano / Enhanced Volcano Plot Dashboard")

# ----------------- 메인 탭: Volcano / EnhancedVolcano -----------------
main_tabs = st.tabs(["🌋 Volcano Plot", "💥 EnhancedVolcano Plot"])
volcano_tab, enhanced_tab = main_tabs

# ----------------- Volcano Plot -----------------
with volcano_tab:
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # Configure
    with configure_tab:
        fc_cutoff_volcano = st.number_input("Fold Change Cutoff (log2)", value=1.0, key="volcano_fc")
        pval_cutoff_volcano = st.number_input("P-value Cutoff", value=0.05, format="%.4f", key="volcano_pval")

    excel_path = "/data/example_data.csv"
    output_svg_volcano = excel_path.replace(".csv", "_volcano.svg")

    # Run
    with run_tab:
        if st.button("Run Volcano Plot"):
            r_code = f"""
library(readr)
library(ggplot2)
data <- read_csv('{excel_path}')
data <- data[!is.na(data$foldchange) & data$foldchange > 0, ]
data$log2FC <- log2(data$foldchange)
fc_cutoff <- {fc_cutoff_volcano}
pval_cutoff <- {pval_cutoff_volcano}
p_cutoff <- -log10(pval_cutoff)
data$group <- "NS"
data$group[data$log2FC > fc_cutoff & -log10(data$pvalue) > p_cutoff] <- "Up"
data$group[data$log2FC < -fc_cutoff & -log10(data$pvalue) > p_cutoff] <- "Down"
color_palette <- c("NS"="gray","Up"="red","Down"="blue")
volcano_plot <- ggplot(data, aes(x=log2FC, y=-log10(pvalue), color=group)) +
    geom_point(size=2, alpha=0.8) +
    scale_color_manual(values=color_palette) +
    geom_vline(xintercept=c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
    geom_hline(yintercept=p_cutoff, linetype="dashed", color="gray") +
    labs(title="Volcano Plot", x="log2 Fold Change", y="-log10 pvalue") +
    theme_minimal()
ggsave(filename='{output_svg_volcano}', plot=volcano_plot, width=8, height=6, dpi=300, device='svg')
"""
            if run_r_script(r_code, output_svg_volcano):
                st.success("Volcano plot generated!")

    # Result
    with result_tab:
        if os.path.exists(output_svg_volcano):
            st.image(output_svg_volcano, caption="Volcano Plot", use_container_width=True)

    # Download
    with download_tab:
        if os.path.exists(output_svg_volcano):
            with open(output_svg_volcano, "rb") as f:
                st.download_button(
                    label="Download Volcano SVG",
                    data=f,
                    file_name="volcano_plot.svg",
                    mime="image/svg+xml"
                )

# ----------------- EnhancedVolcano Plot -----------------
with enhanced_tab:
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # Configure
    with configure_tab:
        fc_cutoff_enhanced = st.number_input("Fold Change Cutoff (log2)", value=1.0, key="enhanced_fc")
        pval_cutoff_enhanced = st.number_input("P-value Cutoff", value=0.05, format="%.4f", key="enhanced_pval")

    output_svg_enhanced = excel_path.replace(".csv", "_enhanced_volcano.svg")

    # Run
    with run_tab:
        if st.button("Run EnhancedVolcano Plot"):
            r_code = f"""
library(readr)
library(EnhancedVolcano)
data <- read_csv('{excel_path}')
data <- data[!is.na(data$foldchange) & data$foldchange > 0, ]
data$log2FC <- log2(data$foldchange)
res <- data.frame(log2FoldChange=data$log2FC, pvalue=data$pvalue)
rownames(res) <- data$Gene_Symbol
res <- res[!is.na(res$log2FoldChange) & is.finite(res$log2FoldChange), ]
svg('{output_svg_enhanced}', width=10, height=8)
EnhancedVolcano(res,
    lab=NA,
    x='log2FoldChange',
    y='pvalue',
    pCutoff={pval_cutoff_enhanced},
    FCcutoff={fc_cutoff_enhanced},
    title='Enhanced Plot',
    subtitle='Sample',
    pointSize=2.5,
    labSize=2.5)
dev.off()
"""
            if run_r_script(r_code, output_svg_enhanced):
                st.success("EnhancedVolcano plot generated!")

    # Result
    with result_tab:
        if os.path.exists(output_svg_enhanced):
            st.image(output_svg_enhanced, caption="EnhancedVolcano Plot", use_container_width=True)

    # Download
    with download_tab:
        if os.path.exists(output_svg_enhanced):
            with open(output_svg_enhanced, "rb") as f:
                st.download_button(
                    label="Download EnhancedVolcano SVG",
                    data=f,
                    file_name="enhanced_volcano_plot.svg",
                    mime="image/svg+xml"
                )
