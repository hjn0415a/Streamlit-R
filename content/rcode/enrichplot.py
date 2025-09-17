import os
import subprocess
import tempfile
import streamlit as st

from src.common.common import page_setup

params = page_setup()

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

a
    """)

result = subprocess.run(
    ["Rscript", r_script_path],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)
st.text(result.stdout)
st.text(result.stderr)

st.success("GO enrichment analysis completed!")
