# hash:sha256:46328663e403f83d4e7983b435d62a234eb044b1615a0a7174157bcb33feebb9
FROM registry.codeocean.com/codeocean/r-studio:r-4.5.1-1.2204.0-rstudio-server-2025.05.1-513-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev=7.81.0-1ubuntu1.23 \
        libfontconfig1-dev=2.13.1-4.2ubuntu5 \
        libfreetype6-dev=2.11.1+dfsg-1ubuntu0.3 \
        libfribidi-dev=1.0.8-2ubuntu3.1 \
        libglpk-dev=5.0-1 \
        libgmp3-dev=2:6.2.1+dfsg-3ubuntu1 \
        libgsl-dev=2.7.1+dfsg-3 \
        libharfbuzz-dev=2.7.4-1ubuntu3.2 \
        libjpeg-dev=8c-2ubuntu10 \
        libpng-dev=1.6.37-3ubuntu0.4 \
        libssl-dev=3.0.2-0ubuntu1.23 \
        libtiff5-dev=4.3.0-6ubuntu0.13 \
        libxml2-dev=2.9.13+dfsg-1ubuntu0.11 \
    && rm -rf /var/lib/apt/lists/*

COPY postInstall /
RUN /postInstall
