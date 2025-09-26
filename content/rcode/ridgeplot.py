import os
import subprocess
import tempfile
import shutil
import streamlit as st

from src.common.common import page_setup
params = page_setup()

st.title("GSEA Ridgeplot")

# ----------------- Main Tabs -----------------
main_tabs = st.tabs(["📊 Ridgeplot (GSEA)"])
ridge_tab = main_tabs[0]

with ridge_tab:

    # ----------------- Sub Tabs -----------------
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Configure -----------------
    with configure_tab:
        input_file = "/data/example_data.csv"  # CSV 입력 경로 (고정)
        output_dir = "/data/ridgeplots"        # 결과 저장 경로 (고정)
        os.makedirs(output_dir, exist_ok=True)

        width = st.number_input("Plot width", value=10.0, step=0.5)
        height = st.number_input("Plot height", value=8.0, step=0.5)

    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run Ridgeplot GSEA"):
            with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
                r_script_path = tmp_r.name
                tmp_r.write(f"""
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

## 1) 경로
file_path <- "{input_file}"
out_dir   <- "{output_dir}"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 2) 데이터
df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("Geneid" %in% names(df))

## 3) 그룹 컬럼 자동 인식
num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
cand <- setdiff(num_cols, c("pvalue","padj","FDR","qvalue","P.Value","p_val","p_val_adj",
                            "log2FC","foldchang","foldchange","foldchge","stat","t","t_stat"))
pattA <- "(^|[^A-Za-z0-9])(A|GroupA|ctrl|control|con|vehicle|veh|untreat|baseline|wt|healthy|pre)($|[^A-Za-z0-9])"
pattB <- "(^|[^A-Za-z0-9])(B|GroupB|case|treated|tx|ko|mut|disease|stim|post|drug)($|[^A-Za-z0-9])"
a_cols <- cand[grepl(pattA, cand, ignore.case = TRUE, perl = TRUE)]
b_cols <- cand[grepl(pattB, cand, ignore.case = TRUE, perl = TRUE)]
stopifnot(length(a_cols) > 1, length(b_cols) > 1)

## 4) 표현행렬/디자인
expr_mat <- as.matrix(df[, c(a_cols, b_cols)])
rownames(expr_mat) <- df$Geneid
if (max(expr_mat, na.rm = TRUE) > 50) expr_mat <- log2(expr_mat + 1)

group  <- factor(c(rep("A", length(a_cols)), rep("B", length(b_cols))))
design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)

## 5) limma moderated t
fit  <- lmFit(expr_mat, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(BvsA = B - A, levels = design)))
tval <- fit2$t[, "BvsA"]
geneList <- sort(tval[is.finite(tval)], decreasing = TRUE)

## rank list 저장
rank_list <- list(
  method   = "limma_moderated_t (BvsA)",
  a_cols   = a_cols,
  b_cols   = b_cols,
  geneList = geneList
)
saveRDS(rank_list, file = file.path(out_dir, "rank_list.rds"))
saveRDS(geneList,  file = file.path(out_dir, "geneList_t.rds"))

## 6) BP/CC/MF 각각 GSEA → ridgeplot + gse 객체 저장
for (ont in c("BP","CC","MF")) {{
  gse <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = ont, minGSSize = 10, maxGSSize = 500,
               pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
  
  saveRDS(gse, file = file.path(out_dir, paste0("gse_", ont, ".rds")))
  
  if (is.null(gse) || nrow(gse@result) == 0) next
  
  p <- ridgeplot(gse, showCategory = 20, fill = "p.adjust", label_format = 40) +
       labs(title = paste("GSEA Ridgeplot (GO:", ont, ")"),
            x = "enrichment distribution") +
       theme_bw()
  
  ggsave(file.path(out_dir, paste0("ridgeplot_", ont, ".svg")),
         p, width = {width}, height = {height}, device = "svg")
}}
""")

            result = subprocess.run(
                ["Rscript", r_script_path],
                capture_output=True,
                text=True,
                encoding="utf-8"
            )
            if result.returncode == 0:
                st.success("Ridgeplot GSEA completed!")
                if result.stdout:
                    st.text(result.stdout)
            else:
                st.error("R script execution failed.")
                if result.stderr:
                    st.text(result.stderr)

    # ----------------- Result -----------------
    with result_tab:
        if os.path.exists(output_dir):
            # Ontology별 탭 생성
            ontology_tabs = st.tabs(["BP", "CC", "MF"])
            for ont_tab, ont in zip(ontology_tabs, ["BP", "CC", "MF"]): 
                with ont_tab:
                    st.subheader(f"Ontology: {ont}")

                    # 해당 ont 이름을 포함한 SVG 파일만 필터링
                    ont_svgs = [f for f in os.listdir(output_dir) if f.endswith(".svg") and ont in f]

                    if ont_svgs:
                        # 2개씩 가로 배치
                        for i in range(0, len(ont_svgs), 2):
                            svg_pair = ont_svgs[i:i+2]
                            cols = st.columns(len(svg_pair))

                            for col, svg_file in zip(cols, svg_pair):
                                with col:
                                    st.markdown(f"**{svg_file}**")
                                    st.image(os.path.join(output_dir, svg_file), use_container_width=True)
                    else:
                        st.info(f"No {ont} ridgeplots found.")
        else:
            st.info("Output directory does not exist.")

    # ----------------- Download -----------------
    with download_tab:
        if os.path.exists(output_dir) and os.listdir(output_dir):
            with tempfile.TemporaryDirectory() as tmpdir:
                zip_path = shutil.make_archive(os.path.join(tmpdir, "ridgeplot_results"), "zip", output_dir)
                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download Ridgeplot Results (ZIP)",
                        data=f,
                        file_name="ridgeplot_results.zip",
                        mime="application/zip"
                    )
        else:
            st.info("No files to download.")
