# Base Ubuntu + Bioconductor image
FROM bioconductor/bioconductor_docker:RELEASE_3_21

ARG DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git build-essential ca-certificates \
    software-properties-common gnupg2 jq cron openssh-client \
    python3 python3-venv python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* 

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install R packages via Bioconductor
RUN Rscript -e "if(!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cran.r-project.org')" && \
    Rscript -e "BiocManager::install(c('EnhancedVolcano', 'clusterProfiler', 'org.Mm.eg.db', 'org.Hs.eg.db', 'enrichplot', 'pheatmap', 'limma', 'pathview', 'RCy3'), update=TRUE, ask=FALSE, dependencies=TRUE)" && \
    Rscript -e "install.packages(c('svglite', 'ggplot2', 'readxl', 'factoextra', 'ggrepel', 'readr', 'cowplot', 'dplyr'), repos='https://cran.r-project.org')"

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
WORKDIR /app
COPY assets/ /app/assets/
COPY content/rcode /app/content/rcode
COPY gdpr_consent/ /app/gdpr_consent
COPY src/ /app/src
COPY app.py /app/app.py
COPY settings.json /app/settings.json
COPY default-parameters.json /app/default-parameters.json
COPY .streamlit/config.toml /app/.streamlit/config.toml

# Add entrypoint script
RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
    echo 'cron' >> /app/entrypoint.sh && \
    echo 'exec /opt/venv/bin/streamlit run /app/app.py --server.port=${PORT:-8501} --server.address=0.0.0.0' >> /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

EXPOSE 8501
ENV PORT=8501

ENTRYPOINT ["/app/entrypoint.sh"]
