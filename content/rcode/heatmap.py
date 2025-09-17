import os
import subprocess
import tempfile
import streamlit as st

from src.common.common import page_setup

params = page_setup()

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
