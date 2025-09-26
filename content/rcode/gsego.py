import os
import subprocess
import tempfile
import shutil
import streamlit as st
import pandas as pd

from src.common.common import page_setup
params = page_setup()

st.title("GSEA Analysis")

# ----------------- Main Tabs -----------------
main_tabs = st.tabs(["🧬 GSEA"])
with main_tabs[0]:

    # ----------------- Sub Tabs -----------------
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Fixed Paths -----------------
    file_path = "/data/example_data.csv"
    out_dir = "/data/GSEA"
    os.makedirs(out_dir, exist_ok=True)

    # ----------------- Configure -----------------
    with configure_tab:
        orgdb = st.selectbox("OrgDb", options=["org.Hs.eg.db", "org.Mm.eg.db"], index=0)
        min_gs_size = st.number_input("Minimum gene set size (minGSSize)", value=10, step=1)
        max_gs_size = st.number_input("Maximum gene set size (maxGSSize)", value=500, step=10)
        pvalue_cutoff = st.number_input("P-value cutoff", value=0.09, step=0.01, format="%.2f")

    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run GSEA Pipeline"):
            with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False, encoding="utf-8") as tmp_r:
                r_script_path = tmp_r.name
                tmp_r.write(f"""
library(limma)
library(clusterProfiler)
library({orgdb})

file_path <- "{file_path}"
out_dir   <- "{out_dir}"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("Geneid" %in% names(df))

num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
cand <- setdiff(num_cols, c("pvalue","padj","FDR","qvalue","P.Value","p_val","p_val_adj",
                            "log2FC","foldchang","foldchange","foldchge","stat","t","t_stat"))
pattA <- "(^|[^A-Za-z0-9])(A|GroupA|ctrl|control|con|vehicle|veh|untreat|baseline|wt|healthy|pre)($|[^A-Za-z0-9])"
pattB <- "(^|[^A-Za-z0-9])(B|GroupB|case|treated|tx|ko|mut|disease|stim|post|drug)($|[^A-Za-z0-9])"
a_cols <- cand[grepl(pattA, cand, ignore.case = TRUE, perl = TRUE)]
b_cols <- cand[grepl(pattB, cand, ignore.case = TRUE, perl = TRUE)]
stopifnot(length(a_cols) > 1, length(b_cols) > 1)

expr_mat <- as.matrix(df[, c(a_cols, b_cols)])
rownames(expr_mat) <- df$Geneid
if (max(expr_mat, na.rm = TRUE) > 50) expr_mat <- log2(expr_mat + 1)

group  <- factor(c(rep("A", length(a_cols)), rep("B", length(b_cols))))
design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)

fit  <- lmFit(expr_mat, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(BvsA = B - A, levels = design)))
tval <- fit2$t[, "BvsA"]
geneList <- sort(tval[is.finite(tval)], decreasing = TRUE)

rank_list <- list(method   = "limma_moderated_t (BvsA)",
                  a_cols   = a_cols,
                  b_cols   = b_cols,
                  geneList = geneList)

saveRDS(rank_list, file = file.path(out_dir, "rank_list.rds"))
saveRDS(geneList,  file = file.path(out_dir, "geneList_t.rds"))

geneList_df <- data.frame(Gene = names(geneList), Score = as.numeric(geneList))
write.csv(geneList_df, file = file.path(out_dir, "geneList_t.csv"), row.names = FALSE)

rank_list_df <- data.frame(Gene = names(rank_list$geneList), Score = as.numeric(rank_list$geneList))
write.csv(rank_list_df, file = file.path(out_dir, "rank_list.csv"), row.names = FALSE)

for (ont in c("BP","CC","MF")) {{
  gse <- gseGO(geneList = geneList, OrgDb = {orgdb}, keyType = "SYMBOL",
               ont = ont, minGSSize = {min_gs_size}, maxGSSize = {max_gs_size},
               pvalueCutoff = {pvalue_cutoff}, pAdjustMethod = "BH", verbose = FALSE)

  saveRDS(gse, file = file.path(out_dir, paste0("gse_", ont, ".rds")))

  if (is.null(gse) || is.null(gse@result) || nrow(gse@result) == 0) {{
    write.csv(data.frame(), file = file.path(out_dir, paste0("gse_", ont, ".csv")), row.names = FALSE)
    next
  }}

  gse_df <- as.data.frame(gse)
  write.csv(gse_df, file = file.path(out_dir, paste0("gse_", ont, ".csv")), row.names = FALSE)
}}
""")

            result = subprocess.run(
                ["Rscript", r_script_path],
                capture_output=True,
                text=True,
                encoding="utf-8"
            )
            if result.returncode == 0:
                st.success("GSEA pipeline completed!")
                if result.stdout:
                    st.text(result.stdout)
            else:
                st.error("R script execution failed.")
                if result.stderr:
                    st.text(result.stderr)

    # ----------------- Result -----------------
    with result_tab:
        ontologies = ["BP", "CC", "MF"]
        ontology_tabs = st.tabs(ontologies)

        for ont_tab, ont in zip(ontology_tabs, ontologies):
            with ont_tab:
                st.markdown(f"**gse_{ont}.csv**")
                csv_file = f"gse_{ont}.csv"
                csv_path = os.path.join(out_dir, csv_file)

                if os.path.exists(csv_path):
                    df = pd.read_csv(csv_path)
                    if df.empty:
                        st.info(f"No results for {ont}.")
                    else:
                        st.dataframe(df)
                else:
                    st.info(f"{csv_file} not found.")

    # ----------------- Download -----------------
    with download_tab:
        if os.listdir(out_dir):
            with tempfile.TemporaryDirectory() as tmpdir:
                zip_path = shutil.make_archive(os.path.join(tmpdir, "GSEA_results"), "zip", out_dir)
                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download GSEA Results (ZIP)",
                        data=f,
                        file_name="GSEA_results.zip",
                        mime="application/zip"
                    )
        else:
            st.info("No files to download.")
