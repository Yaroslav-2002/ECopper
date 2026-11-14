#!/bin/sh
# Prepare the database and BLAST indexes, then start the given command (gunicorn by default).
set -e

# A secret key generated on first start and kept in the volume, unless one is passed in.
if [ -z "${DJANGO_SECRET_KEY}" ]; then
    key_file=/app/var/secret_key
    if [ ! -s "${key_file}" ]; then
        python -c "import secrets; print(secrets.token_urlsafe(50))" > "${key_file}"
    fi
    DJANGO_SECRET_KEY="$(cat "${key_file}")"
    export DJANGO_SECRET_KEY
fi

python manage.py migrate --noinput
# Rebuild the copper reference table and BLAST protein database from data/ (takes a few seconds).
python manage.py build_cu_db

exec "$@"
