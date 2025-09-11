# Base Ubuntu image
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git build-essential ca-certificates \
    software-properties-common gnupg2 jq cron openssh-client && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Miniforge (conda)
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh
ENV PATH="/opt/conda/bin:$PATH"

# Create Conda environment
RUN conda install -y mamba -n base -c conda-forge && \
    mamba create -y -n quantms-env python=3.10 && \
    mamba install -y -n quantms-env -c conda-forge \
        r-base rpy2 r-svglite r-ggplot2 r-readxl && \
    conda clean -afy

RUN conda run -n quantms-env Rscript -e "if(!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cran.r-project.org')" && \
    conda run -n quantms-env Rscript -e "BiocManager::install(c('EnhancedVolcano', 'clusterProfiler', 'org.Mm.eg.db', 'enrichplot', 'pheatmap'), update=FALSE, ask=FALSE)" && \
    # conda run -n quantms-env Rscript -e "install.packages(c('factoextra', 'ggrepel'), repos='https://cran.r-project.org')"

# Install Python dependencies
COPY requirements.txt .
SHELL ["conda", "run", "-n", "quantms-env", "/bin/bash", "-c"]
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
WORKDIR /app
COPY assets/ /app/assets/
COPY content/quantms /app/content/quantms
COPY gdpr_consent/ /app/gdpr_consent
COPY src/ /app/src
COPY app.py /app/app.py
COPY settings.json /app/settings.json
COPY default-parameters.json /app/default-parameters.json
COPY .streamlit/config.toml /app/.streamlit/config.toml

# Add entrypoint script
SHELL ["/bin/bash", "-c"]
RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
   echo 'cron' >> /app/entrypoint.sh && \
   echo 'exec conda run -n quantms-env streamlit run /app/app.py --server.port=${PORT:-8501} --server.address=0.0.0.0' >> /app/entrypoint.sh && \
   chmod +x /app/entrypoint.sh

EXPOSE 8501
ENV PORT=8501

ENTRYPOINT ["/app/entrypoint.sh"]