FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
ENV HOME=/home

# Set default entrypoint for the build to be the bash shell
ENTRYPOINT ["bash"]

# Install a few ubuntu dependencies
RUN apt-get update && \
    apt install --fix-broken && \
    apt-get install -y \
    build-essential \
    curl \
    wget \
    make \
    gcc \
    cmake \
    libxml2-dev \
    libxslt-dev \
    libffi-dev \
    git && \
    apt install --fix-broken && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /dependencies && \
    dpkg -l > /dependencies/apt-get.lock

# Install everything else with Pixi:
# --------------------------------
# 1) copy the required dependency and configuration file into the image
COPY pyproject.toml $HOME/pyproject.toml
COPY pixi.lock $HOME/pixi.lock

# 2) install pixi
RUN cd $HOME && PIXI_ARCH=x86_64 curl -fsSL https://pixi.sh/install.sh | bash

# 3) make sure pixi and pixi installs are on the $PATH
ENV PATH=$PATH:$HOME/.pixi/bin

# 4) install everything else with pixi
RUN cd $HOME && pixi install --verbose --color=always --frozen && pixi clean cache --assume-yes

# 5) modify the shell config so that each container launches within the pixi env
RUN echo "export PATH=$PATH:$HOME/.pixi/envs/default/bin" >> $HOME/.bashrc

# 6) modify some nextflow environment variables
RUN echo "export NXF_CACHE_DIR=/scratch" >> $HOME/.bashrc
RUN echo "export NXF_HOME=/scratch" >> $HOME/.bashrc

