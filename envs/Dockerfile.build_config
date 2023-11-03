FROM python:3.10

COPY ./requirements-build_config.txt /requirements.txt
RUN apt-get update && \
    apt-get install -y --no-install-recommends && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install --upgrade pip --no-cache-dir && \
    pip install --no-cache-dir -r /requirements.txt
