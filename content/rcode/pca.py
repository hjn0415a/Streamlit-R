import os
import subprocess
import tempfile
import streamlit as st

from src.common.common import page_setup

params = page_setup()

st.markdown("## PCA Plot (factoextra + ggrepel)")

excel_path_pca = "/data/proteinGroups.xlsx"  # 컨테이너 안 경로
output_svg_pca = "/data/PCA_plot6.svg"

with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
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
