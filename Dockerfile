################## Base Image ##########
FROM --platform=linux/amd64 debian:latest 

################## ARGUMENTS/Environments ##########
ARG BUILD_DATE
ARG BUILD_VERSION
ARG LICENSE="Apache-2.0"
ARG PYTHON_VERSION="3.6.10"
ARG POSTPROCESSING_VARIANT_CALLS_VERSION="re-project-layout"
ARG VCS_REF
ENV BOOST_ROOT /usr
ENV PATH="/root/miniconda3/bin:$PATH"
ARG PATH="/root/miniconda3/bin:$PATH"

################## METADATA ########################
LABEL org.opencontainers.image.vendor="MSKCC"
LABEL org.opencontainers.image.authors="Eric Buehlere (buehlere@mskcc.org)"

LABEL org.opencontainers.image.created=${BUILD_DATE} \
    org.opencontainers.image.version=${BUILD_VERSION} \
    org.opencontainers.image.licenses=${LICENSE} \
    org.opencontainers.image.version.python=${PYTHON_VERSION} \ 
    org.opencontainers.image.vcs-ref=${VCS_REF}

LABEL org.opencontainers.image.description="This container uses conda/conda/miniconda3 as the base image to build"

################## INSTALL ##########################
# installing tools 
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget 

# install postprocessing_variant_calls 
RUN cd /opt \ 
    && git clone --recursive -b ${POSTPROCESSING_VARIANT_CALLS_VERSION} https://github.com/msk-access/postprocessing_variant_calls.git \ 
    && cd postprocessing_variant_calls

# install miniconda and add to path 
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-4.5.4-Linux-x86_64.sh -b \
    && rm -f Miniconda3-4.5.4-Linux-x86_64.sh 

RUN conda install -c bioconda pyvcf

RUN cd /opt/postprocessing_variant_calls \
    && make deps-install