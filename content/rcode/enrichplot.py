import os
import subprocess
import tempfile
import shutil
import streamlit as st
import pandas as pd

from src.common.common import page_setup
params = page_setup()

st.title("Enrichplot")

# ----------------- GO Enrichment -----------------
main_tabs = st.tabs(["🧬 GO Enrichment"])
enrich_tab = main_tabs[0]

with enrich_tab:
    sub_tabs = st.tabs(["⚙️ Configure", "🚀 Run", "📊 Result", "⬇️ Download"])
    configure_tab, run_tab, result_tab, download_tab = sub_tabs

    # ----------------- Configure -----------------
    with configure_tab:
        result_root = "/data/Deg"                 # DEG 결과 폴더
        output_root = "/data/Enrichment/out"      # Enrichment 결과 저장 폴더
        os.makedirs(output_root, exist_ok=True)

        showCategory = st.number_input(
            "Number of categories to show (showCategory)", value=10, step=1
        )
        pvalueCutoff = st.number_input(
            "P-value cutoff", value=0.9, step=0.01, format="%.3f"
        )
        org_db = st.selectbox(
            "OrgDb", options=["org.Hs.eg.db", "org.Mm.eg.db"], index=0
        )
        plot_width = st.number_input("Plot width", value=8.0, step=0.5)
        plot_height = st.number_input("Plot height", value=6.0, step=0.5)

    # ----------------- Run -----------------
    with run_tab:
        if st.button("Run GO Enrichment"):
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".R", delete=False, encoding="utf-8"
            ) as tmp_r:
                r_script_path = tmp_r.name
                tmp_r.write(f"""
library(clusterProfiler)
library({org_db})
library(enrichplot)
library(ggplot2)

run_enrich_genedi_min <- function(result_root,
                                  output_root,
                                  combo_names,
                                  file_name    = "filtered_gene_list.csv",
                                  showCategory = {showCategory},
                                  p_cut        = {pvalueCutoff},
                                  save_ego     = TRUE,
                                  width        = {plot_width},
                                  height       = {plot_height}) {{

    for (nm in combo_names) {{
        combo_dir_in <- file.path(result_root, nm)
        f <- file.path(combo_dir_in, file_name)
        if (!file.exists(f)) next

        df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)

        sym_col <- grep("^(Geneid|Gene_Symbol|SYMBOL)$", names(df),
                        ignore.case = TRUE, value = TRUE)[1]
        if (is.na(sym_col)) next

        conv <- tryCatch(
            bitr(df[[sym_col]], fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = {org_db}),
            error = function(e) {{ NULL }}
        )
        if (is.null(conv) || !"ENTREZID" %in% names(conv)) next

        ids <- unique(na.omit(conv$ENTREZID))
        if (!length(ids)) next

        combo_dir_out <- file.path(output_root, nm)
        fig_dir <- file.path(combo_dir_out, "figure")
        if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

        for (ont in c("BP", "CC", "MF")) {{
            ego <- suppressMessages(
                enrichGO(
                    gene           = ids,
                    OrgDb          = {org_db},
                    keyType        = "ENTREZID",
                    ont            = ont,
                    pAdjustMethod  = "BH",
                    pvalueCutoff   = p_cut,
                    qvalueCutoff   = 1,
                    readable       = TRUE
                )
            )

            if (!is.null(ego) && !is.null(ego@result) && nrow(ego@result) > 0) {{
                # 테이블 저장
                write.csv(ego@result,
                          file.path(combo_dir_out, sprintf("GO_%s_result.csv", ont)),
                          row.names = FALSE)

                # 플롯 저장
                p <- dotplot(ego, showCategory = showCategory,
                             x = "GeneRatio", color = "p.adjust") +
                     ggtitle(sprintf("GO %s - %s", ont, nm))
                ggsave(file.path(fig_dir, sprintf("GO_%s.svg", ont)),
                       p, width = width, height = height)

                # ego 객체 저장
                if (isTRUE(save_ego)) {{
                    saveRDS(ego, file.path(combo_dir_out, sprintf("GO_%s_ego.rds", ont)))
                }}
            }}
        }}
    }}
}}

combo_names <- readRDS(file.path("{result_root}", "combo_names.rds"))

run_enrich_genedi_min(
    result_root = "{result_root}",
    output_root = "{output_root}",
    combo_names = combo_names,
    file_name   = "filtered_gene_list.csv",
    showCategory= {showCategory},
    p_cut       = {pvalueCutoff},
    save_ego    = TRUE,
    width       = {plot_width},
    height      = {plot_height}
)
""")

            result = subprocess.run(
                ["Rscript", r_script_path],
                capture_output=True,
                text=True,
                encoding="utf-8"
            )
            if result.returncode == 0:
                st.success("GO enrichment analysis completed!")
                st.text(result.stdout)
            else:
                st.error("R script execution failed.")
                st.text(result.stderr)

    # ----------------- Result -----------------
    with result_tab:
        combo_csv = os.path.join(result_root, "combo_names.csv")
        if os.path.exists(combo_csv):
            combos = pd.read_csv(combo_csv)["combo"].tolist()

            # Ontology 별 탭 생성
            ontology_tabs = st.tabs(["BP", "CC", "MF"])
            for ont_tab, ont in zip(ontology_tabs, ["BP", "CC", "MF"]):
                with ont_tab:
                    st.subheader(f"Ontology: {ont}")

                    # 각 조합별 결과 처리
                    for i in range(0, len(combos), 2):
                        combo_pair = combos[i:i+2]  # 최대 2개씩 열 배치
                        cols = st.columns(len(combo_pair))
    
                        for col, combo in zip(cols, combo_pair):
                            with col:
                                # 조합 이름과 유전자 수 표시
                                result_file = os.path.join(output_root, combo, f"GO_{ont}_result.csv")
                                plot_file   = os.path.join(output_root, combo, "figure", f"GO_{ont}.svg")

                                if os.path.exists(result_file):
                                    df = pd.read_csv(result_file)
                                    st.markdown(f"**{combo}** (Genes: {len(df)})")
                                
                                    # 이미지 먼저 출력
                                    if os.path.exists(plot_file):
                                        st.image(plot_file, use_container_width=True)
                                
                                    # 테이블 출력
                                    st.dataframe(df, use_container_width=True, height=250)
                                else:
                                    st.markdown(f"**{combo}** - No results for {ont}")

    # ----------------- Download -----------------
    with download_tab:
        combo_csv = os.path.join(result_root, "combo_names.csv")
        if os.path.exists(combo_csv):
            combos = pd.read_csv(combo_csv)["combo"].tolist()
            if combos:
                with tempfile.TemporaryDirectory() as tmpdir:
                    for combo in combos:
                        src = os.path.join(output_root, combo)
                        dst = os.path.join(tmpdir, combo)
                        if os.path.exists(src):
                            shutil.copytree(src, dst)
                    zip_path = shutil.make_archive(os.path.join(tmpdir, "Enrichment_combos"), "zip", tmpdir)
                    with open(zip_path, "rb") as f:
                        st.download_button(
                            label="Download Enrichment Results (ZIP)",
                            data=f,
                            file_name="Enrichment_combos.zip",
                            mime="application/zip"
                        )
