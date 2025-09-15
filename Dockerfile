# Base Ubuntu + Bioconductor image
FROM bioconductor/bioconductor_docker

ARG DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git build-essential ca-certificates \
    software-properties-common gnupg2 jq cron openssh-client \
    python3 python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages via Bioconductor
RUN Rscript -e "if(!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cran.r-project.org')" && \
    Rscript -e "BiocManager::install(c('EnhancedVolcano', 'clusterProfiler', 'org.Mm.eg.db', 'enrichplot', 'pheatmap'), update=TRUE, ask=FALSE, dependencies=TRUE)" && \
    Rscript -e "install.packages(c('rpy2', 'svglite', 'ggplot2', 'readxl'), repos='https://cran.r-project.org')"

# Install Python dependencies
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

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
RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
    echo 'cron' >> /app/entrypoint.sh && \
    echo 'exec streamlit run /app/app.py --server.port=${PORT:-8501} --server.address=0.0.0.0' >> /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

EXPOSE 8501
ENV PORT=8501

ENTRYPOINT ["/app/entrypoint.sh"]