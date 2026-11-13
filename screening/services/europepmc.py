import requests
from django.conf import settings

from ..models import LiteratureCache


def hit_count(gene_name):
    key = gene_name.strip()
    if not key:
        return None
    cached = LiteratureCache.objects.filter(gene_name=key).first()
    if cached:
        return cached.hit_count
    params = {"query": f'"{key}" AND bioremediation', "format": "json", "pageSize": 1}
    try:
        resp = requests.get(settings.EUROPEPMC_URL, params=params, timeout=settings.EUROPEPMC_TIMEOUT)
        resp.raise_for_status()
        count = int(resp.json().get("hitCount", 0))
    except (requests.RequestException, ValueError, TypeError):
        return None
    LiteratureCache.objects.update_or_create(gene_name=key, defaults={"hit_count": count})
    return count
