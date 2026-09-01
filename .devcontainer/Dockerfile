FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        r-base \
        r-base-dev \
        gfortran \
        libcurl4-openssl-dev \
        libgit2-dev \
        libssl-dev \
        libxml2-dev \
        make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*