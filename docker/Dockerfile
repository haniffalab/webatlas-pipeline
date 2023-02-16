FROM python:3.9
ARG BIOFORMATS2RAW="0.4.0"

COPY ./requirements.txt /requirements.txt
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      wget unzip cmake g++ ant procps libblosc1 libvips libtiff-tools openjdk-17-jre-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    wget -q -O bioformats2raw.zip https://github.com/glencoesoftware/bioformats2raw/releases/download/v${BIOFORMATS2RAW}/bioformats2raw-${BIOFORMATS2RAW}.zip && \
    unzip bioformats2raw.zip -d /usr/local/share/ && \
    ln -sf /usr/local/share/bioformats2raw-${BIOFORMATS2RAW}/bin/bioformats2raw /usr/local/bin/bioformats2raw && \
    rm -rf bioformats2raw.zip

RUN pip install --upgrade pip setuptools distlib --no-cache-dir && \
    pip install --no-cache-dir -r /requirements.txt

ENTRYPOINT []