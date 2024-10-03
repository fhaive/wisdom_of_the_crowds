FROM ubuntu:24.04

RUN apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
        wget bzip2 ca-certificates \
        git \
        tini \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget --no-hsts "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -O /opt/miniforge.sh \
    && /bin/bash /opt/miniforge.sh -b -p /opt/conda \
    && rm /opt/miniforge.sh
#    && echo ". /opt/conda/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc \
#    && echo ". /opt/conda/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

ENV PATH="/opt/conda/bin:$PATH"

#SHELL ["/bin/bash", "--login"]
RUN conda init bash \
    && conda config --set solver libmamba \
    && conda create -n renv -y -c conda-forge -c r -c bioconda r-base=4.2.2 r-data.table=1.14.10 r-dt=0.27 bioconductor-complexheatmap=2.14.0 r-circlize=0.4.15 bioconductor-interactivecomplexheatmap=1.6.0 r-kableextra=1.3.4 r-plotly=4.10.1 r-dplyr=1.1.3 r-ggplot2=3.4.2 r-reshape2=1.4.4 r-magrittr=2.0.3 r-readxl=1.4.2 r-shiny=1.7.4 r-shinybs=0.61.1 r-shinyjs=2.1.0 r-shinywidgets=0.7.6 r-shinydashboard=0.7.2

COPY app/ /app/
COPY runapp.R /runapp.R

CMD conda run -n renv --no-capture-output Rscript /runapp.R
