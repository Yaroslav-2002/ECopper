# CuRemediate: Django app + NCBI BLAST+ in one image.
# Build and run with:  docker compose up --build

# --- Stage 1: download the official BLAST+ release and keep only the three programs the app uses.
FROM debian:bookworm-slim AS blast
ARG BLAST_VERSION=2.17.0
ARG TARGETARCH
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN set -eux; \
    case "${TARGETARCH:-amd64}" in \
      amd64) arch=x64 ;; \
      arm64) arch=aarch64 ;; \
      *) echo "Unsupported architecture: ${TARGETARCH}" >&2; exit 1 ;; \
    esac; \
    file="ncbi-blast-${BLAST_VERSION}+-${arch}-linux.tar.gz"; \
    url="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/${file}"; \
    cd /tmp; \
    curl -fsSLO "${url}"; \
    curl -fsSL "${url}.md5" | md5sum -c -; \
    tar -xzf "${file}"; \
    mkdir -p /opt/blast/bin; \
    for prog in blastp blastx makeblastdb; do \
      cp "ncbi-blast-${BLAST_VERSION}+/bin/${prog}" /opt/blast/bin/; \
    done; \
    rm -rf /tmp/*

# --- Stage 2: the application.
FROM python:3.11-slim-bookworm
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    BLAST_BIN_DIR=/opt/blast/bin \
    CUREMEDIATE_DB_PATH=/app/var/db.sqlite3 \
    DJANGO_DEBUG=0
# libgomp1: OpenMP runtime used by the BLAST+ programs.
RUN apt-get update \
 && apt-get install -y --no-install-recommends libgomp1 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=blast /opt/blast /opt/blast

WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt "gunicorn>=22"

COPY . .
RUN sed -i 's/\r$//' docker/entrypoint.sh \
 && DJANGO_SECRET_KEY=collectstatic-only python manage.py collectstatic --noinput \
 && useradd --create-home --uid 1000 app \
 && mkdir -p /app/var \
 && chown -R app:app /app/var /app/data
USER app

VOLUME /app/var
EXPOSE 8000
ENTRYPOINT ["sh", "/app/docker/entrypoint.sh"]
# BLAST on a large batch can take minutes, hence the long worker timeout.
CMD ["gunicorn", "config.wsgi:application", "--bind", "0.0.0.0:8000", "--workers", "2", "--timeout", "600"]
