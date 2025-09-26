import os
import subprocess
import tempfile
import streamlit as st

from src.common.common import page_setup

params = page_setup()

st.title("Heatmap")

# ----------------- Heatmap -----------------
main_tabs = st.tabs(["🌡️ Heatmap"])
heatmap_tab = main_tabs[0]

with heatmap_tab:
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Configure -----------------
    with configure_tab:
        width_heatmap = st.number_input("Plot Width", value=8.0, step=0.5)
        height_heatmap = st.number_input("Plot Height", value=10.0, step=0.5)
        top_n_genes = st.number_input("Top N genes (by p-value)", value=50, step=5)

    # 입력/출력 경로
    csv_path_heatmap = "/data/example_data.csv"
    output_svg_heatmap = "/data/heatmap2.svg"

    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run Heatmap"):
            r_code = f"""
library(pheatmap)
library(readr)
library(svglite)

# CSV 데이터 불러오기
data <- read_csv("{csv_path_heatmap}")

# 행 이름 설정 (첫 번째 열: geneid)
gene_names <- data[[1]]
data <- data[, -1]

# 샘플 컬럼만 추출 (숫자로 끝나는 열만)
sample_cols <- grep("([0-9]+$)", names(data), value = TRUE)

# 행렬 생성 (샘플 열만 사용)
mat <- as.matrix(data[, sample_cols, drop = FALSE])
rownames(mat) <- gene_names

# 그룹 자동 추출
annotation_col <- data.frame(
  Group = factor(sub("(_[0-9]+$)|([0-9]+$)", "", sample_cols)),
  row.names = sample_cols
)

# SVG 파일로 저장
svglite("{output_svg_heatmap}", width = {width_heatmap}, height = {height_heatmap})

# 히트맵 생성 (p-value 기준 상위 N개 추출)
pheatmap(mat[order(data$pvalue)[1:{top_n_genes}], , drop = FALSE],
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         color = colorRampPalette(c("#6699e0", "white", "#e06666"))(100)
)

dev.off()
"""
            try:
                subprocess.run(["Rscript", "-e", r_code], check=True, text=True)
                st.success("Heatmap generated successfully!")
            except subprocess.CalledProcessError as e:
                st.error(f"Heatmap generation failed: {e}")

    # ----------------- Result -----------------
    with result_tab:
        if os.path.exists(output_svg_heatmap):
            st.image(output_svg_heatmap, caption="Heatmap", use_container_width=True)

    # ----------------- Download -----------------
    with download_tab:
        if os.path.exists(output_svg_heatmap):
            with open(output_svg_heatmap, "rb") as f:
                st.download_button(
                    label="Download Heatmap SVG",
                    data=f,
                    file_name="heatmap.svg",
                    mime="image/svg+xml"
                )
