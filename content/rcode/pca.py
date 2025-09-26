import os
import subprocess
import tempfile
import streamlit as st

from src.common.common import page_setup

params = page_setup()

st.title("PCA (Principal Component Analysis)")

# ----------------- PCA -----------------
main_tabs = st.tabs(["📉 PCA Plot"])
pca_tab = main_tabs[0]

with pca_tab:
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Configure -----------------
    with configure_tab:
        width_pca = st.number_input("Plot Width", value=8.0, step=0.5, key="pca_width")
        height_pca = st.number_input("Plot Height", value=6.0, step=0.5, key="pca_height")
        pointshape_pca = st.number_input("Point Shape", value=16, step=1, key="pca_pointshape")
        pointsize_pca = st.number_input("Point Size", value=3.5, step=0.5, key="pca_pointsize")
        text_size_pca = st.number_input("Label Text Size", value=4.0, step=0.5, key="pca_textsize")

    # 입력/출력 경로
    csv_path_pca = "/data/example_data.csv"
    output_svg_pca = "/data/PCA_plot2.svg"

    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run PCA Plot"):
            with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
                r_script_path = tmp_r.name
                tmp_r.write(f"""
library(readr)
library(factoextra)
library(ggrepel)
library(svglite)

# 1) 데이터 읽기
dat <- as.data.frame(read_csv("{csv_path_pca}"))

# 2) 샘플(그룹) 열만 추출 (Geneid, foldchange, pvalue 제외)
sample_cols <- grep("(^Group)|(_[0-9]+$)", names(dat), value = TRUE)

# 3) (샘플 × 유전자) 행렬 생성
X <- t(as.matrix(dat[, sample_cols, drop = FALSE]))
rownames(X) <- sample_cols

# 4) Z-score 정규화
Xz <- scale(X, center = TRUE, scale = TRUE)

# 5) 결측치(NA) 있는 열 제거
Xz <- Xz[, colSums(is.na(Xz)) == 0, drop = FALSE]

# 6) PCA
pca_res <- prcomp(Xz, center = FALSE, scale. = FALSE)
rownames(pca_res$x) <- rownames(X)

# 7) 그룹 자동 추출
sample_groups <- factor(sub("(_[0-9]+$)|([0-9]+$)", "", rownames(pca_res$x)))

# 8) 데이터프레임 생성
df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                 sample = rownames(pca_res$x), group = sample_groups)

# 9) SVG 파일로 저장
svglite("{output_svg_pca}", width = {width_pca}, height = {height_pca})

p <- fviz_pca_ind(
  pca_res,
  geom.ind   = "point",
  col.ind    = sample_groups,
  pointshape = {pointshape_pca},
  pointsize  = {pointsize_pca},
  mean.point = FALSE,
  addEllipses= FALSE
) +
  geom_text_repel(data = df, aes(PC1, PC2, label = sample, color = group),
                  size = {text_size_pca}, show.legend = FALSE)

print(p)
dev.off()
                """)

            try:
                subprocess.run(["Rscript", r_script_path], check=True, text=True)
                st.success("PCA plot generated successfully!")
            except subprocess.CalledProcessError as e:
                st.error(f"PCA plot generation failed: {e}")
            finally:
                os.remove(r_script_path)

    # ----------------- Result -----------------
    with result_tab:
        if os.path.exists(output_svg_pca):
            st.image(output_svg_pca, caption="PCA Plot", use_container_width=True)

    # ----------------- Download -----------------
    with download_tab:
        if os.path.exists(output_svg_pca):
            with open(output_svg_pca, "rb") as f:
                st.download_button(
                    label="Download PCA SVG",
                    data=f,
                    file_name="PCA_plot.svg",
                    mime="image/svg+xml"
                )
