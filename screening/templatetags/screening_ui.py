"""Display helpers for the screening templates."""
from django import template
from django.utils.html import conditional_escape
from django.utils.safestring import mark_safe

from ..services.scoring import NO_EVIDENCE, NOT_TRANSFERABLE, format_evalue

register = template.Library()

CATEGORY_CLASSES = {
    "High potential": "high",
    "Moderate": "moderate",
    "Low": "low",
    NOT_TRANSFERABLE: "nt",
    NO_EVIDENCE: "none",
}


@register.filter
def evalue(value, digits=1):
    """E-value in scientific notation; BLAST's 0.0 (too small to print) becomes '< 1e-180'."""
    return format_evalue(value, int(digits))


@register.filter
def breakable(text):
    """Allow long sequence IDs (sp|P36649|CUEO_ECOLI) to wrap after '|' and ':' instead of mid-word."""
    escaped = conditional_escape(text)
    return mark_safe(escaped.replace("|", "|<wbr>").replace(":", ":<wbr>"))


@register.filter
def category_class(category):
    """CSS modifier for a ranking category; '' (old scoring scheme) maps to 'legacy'."""
    return CATEGORY_CLASSES.get(category, "legacy")


@register.filter
def uncertainty_class(uncertainty):
    return {"Low": "good", "Moderate": "mid", "High": "bad"}.get(uncertainty, "")


@register.filter
def segments(points, maximum):
    """One boolean per available point, True for each point earned (points None = skipped criterion)."""
    earned = points or 0
    return [i < earned for i in range(int(maximum or 0))]


@register.filter
def percent_of(value, total):
    """value as a whole-number percentage of total (0 when total is empty)."""
    try:
        return round(100 * float(value) / float(total)) if total else 0
    except (TypeError, ValueError):
        return 0
