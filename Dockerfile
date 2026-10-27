# This Dockerfile was crested to complete test task
# Базовый образ Ubuntu:22.04
FROM ubuntu:22.04

# Отключить интерактивные запросы при установке пакетов
ENV DEBIAN_FRONTEND=noninteractive

# Переменная SOFT согласно заданию
ENV SOFT=/soft

# Слой 1: Установка общих (не специализированных) пакетов
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    curl \
    cmake \
    ca-certificates \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libgsl-dev \
    perl \
    cpanminus \
    autoconf \
    automake \
    libtool \
    pkg-config \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Слой 2: Установка pysam.
# pysam. Версия: 0.24.1 (дата релиза: 2026-09-07, согласно данным https://github.com/pysam-developers/pysam/releases?ysclid=mvjydhxrq555073988)
RUN pip3 install --no-cache-dir pysam==0.24.1

# Слой 3: Установка libdeflate. Используется samtools/htslib для ускорения сжатия.
# libdeflate. Версия: 1.26 (дата релиза: 2026-08-22, согласно данным https://github.com/ebiggers/libdeflate/releases)
WORKDIR /tmp
RUN wget https://github.com/ebiggers/libdeflate/archive/refs/tags/v1.26.tar.gz && \
    tar -xzf v1.26.tar.gz && \
    cd libdeflate-1.26 && \
    cmake -B build -DCMAKE_INSTALL_PREFIX=$SOFT/libdeflate-1.26 && \
    cmake --build build -j$(nproc) && \
    cmake --install build && \
    cd /tmp && rm -rf libdeflate-1.26 v1.26.tar.gz

# Добавление пути к libdeflate в переменные окружения для сборки htslib.
ENV CPATH=$SOFT/libdeflate-1.26/include
ENV LIBRARY_PATH=$SOFT/libdeflate-1.26/lib
ENV LD_LIBRARY_PATH=$SOFT/libdeflate-1.26/lib:${LD_LIBRARY_PATH:-}

# Слой 4: Установка HTSlib. Требуется для samtools и bcftools.
# HTSlib. Версия: 1.24 (дата релиза: 2026-07-09, согласно данным https://github.com/samtools/htslib/releases)
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.24/htslib-1.24.tar.bz2 && \
    tar -xjf htslib-1.24.tar.bz2 && \
    cd htslib-1.24 && \
    ./configure --prefix=$SOFT/htslib-1.24 \
                --with-libdeflate \
                --enable-libcurl && \
    make -j$(nproc) && \
    make install && \
    cd /tmp && rm -rf htslib-1.24 htslib-1.24.tar.bz2

# Слой 5: Установка SAMtools.
# SAMtools. Версия: 1.24 (дата релиза: 2026-07-09, согласно данным https://github.com/samtools/samtools/releases)
WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.24/samtools-1.24.tar.bz2 && \
    tar -xjf samtools-1.24.tar.bz2 && \
    cd samtools-1.24 && \
    ./configure --prefix=$SOFT/samtools-1.24 \
                --with-htslib=$SOFT/htslib-1.24 && \
    make -j$(nproc) && \
    make install && \
    cd /tmp && rm -rf samtools-1.24 samtools-1.24.tar.bz2

# Слой 6: Установка BCFtools.
# BCFtools. Версия: 1.24 (дата релиза: 2026-07-09, согласно данным https://github.com/samtools/bcftools/releases)
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/1.24/bcftools-1.24.tar.bz2 && \
    tar -xjf bcftools-1.24.tar.bz2 && \
    cd bcftools-1.24 && \
    ./configure --prefix=$SOFT/bcftools-1.24 \
                --with-htslib=$SOFT/htslib-1.24 && \
    make -j$(nproc) && \
    make install && \
    cd /tmp && rm -rf bcftools-1.24 bcftools-1.24.tar.bz2

# Слой 7: Установка VCFtools.
# VCFtools. Версия: 0.1.17 (дата релиза: 2025-05-15, согласно данным https://github.com/vcftools/vcftools/releases)
WORKDIR /tmp
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.17/vcftools-0.1.17.tar.gz && \
    tar -xzf vcftools-0.1.17.tar.gz && \
    cd vcftools-0.1.17 && \
    ./configure --prefix=$SOFT/vcftools-0.1.17 && \
    make -j$(nproc) && \
    make install && \
    cd /tmp && rm -rf vcftools-0.1.17 vcftools-0.1.17.tar.gz

# Добавление python скрипта
COPY add_ref_allele.py /usr/local/bin/add_ref_allele.py
RUN chmod +x /usr/local/bin/add_ref_allele.py

# Полные пути до исполняемых файлов в PATH
ENV PATH=$SOFT/samtools-1.24/bin:$SOFT/bcftools-1.24/bin:$SOFT/vcftools-0.1.17/bin:$PATH
ENV LD_LIBRARY_PATH=$SOFT/htslib-1.24/lib:$SOFT/libdeflate-1.26/lib:$LD_LIBRARY_PATH

# Отдельные переменные окружения для каждой программы
ENV SAMTOOLS=$SOFT/samtools-1.24/bin/samtools
ENV BCFTOOLS=$SOFT/bcftools-1.24/bin/bcftools
ENV VCFTOOLS=$SOFT/vcftools-0.1.17/bin/vcftools
