"""Sphinx configuration for the public Silex documentation."""

from __future__ import annotations

import os


project = "Silex"
copyright = "2026, William Youmans and contributors"
author = "William Youmans and contributors"

release = os.environ.get("SILEX_DOCS_RELEASE", "0.1.1")
version_parts = release.split(".")
version = ".".join(version_parts[:2]) if len(version_parts) >= 2 else release

docs_channel = os.environ.get("SILEX_DOCS_CHANNEL", "dev").strip().lower()
if docs_channel not in {"stable", "dev"}:
    raise ValueError(
        "SILEX_DOCS_CHANNEL must be either 'stable' or 'dev', "
        f"not {docs_channel!r}"
    )
tags.add(docs_channel)

extensions = []
exclude_patterns = []

rst_epilog = f"""
.. |silex-release| replace:: {release}
.. |docs-channel| replace:: {docs_channel}
"""

html_theme = "alabaster"
html_title = f"Silex {release} ({docs_channel})"
html_short_title = f"Silex {release}"
html_baseurl = (
    "https://wjyoumans.github.io/silex/"
    + (f"{release}/" if docs_channel == "stable" else "dev/")
)
html_theme_options = {
    "description": (
        "Native C++20 computational algebraic number theory over FLINT"
    ),
    "fixed_sidebar": True,
    "github_banner": True,
    "github_button": True,
    "github_count": False,
    "github_repo": "silex",
    "github_type": "star",
    "github_user": "wjyoumans",
}
html_show_sourcelink = True
html_context = {
    "docs_channel": docs_channel,
    "silex_release": release,
}
