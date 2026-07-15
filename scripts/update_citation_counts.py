#!/usr/bin/env python3
"""Refresh the citation-count table in README.md using OpenAlex.

OpenAlex (https://openalex.org) is used instead of Google Scholar because it
has a free, stable, official API. Google Scholar has no API and actively
blocks automated/datacenter-IP requests, which makes it unsuitable for an
unattended, scheduled job like this one. Counts will run a little lower than
Google Scholar's, since Scholar indexes a broader (and less curated) set of
sources.

Uses only the standard library so it can run on a bare GitHub Actions runner
with no extra setup.
"""
import datetime
import json
import pathlib
import urllib.request

README_PATH = pathlib.Path(__file__).resolve().parent.parent / "README.md"
START_MARKER = "<!-- CITATION-STATS:START -->"
END_MARKER = "<!-- CITATION-STATS:END -->"

# (short label, doi) for every paper in the References section, in the same
# order they appear there.
PAPERS = [
    ("MSstats (2014)", "10.1093/bioinformatics/btu305"),
    ("MSstats 4.0 (2023)", "10.1021/acs.jproteome.2c00834"),
    ("MSstats + FragPipe DIA workflow (2024)", "10.1038/s41596-024-01000-3"),
    ("MSstatsTMT (2020)", "10.1074/mcp.RA120.002105"),
    ("MSstatsTMT repeated measures (2023)", "10.1021/acs.jproteome.3c00155"),
    ("MSstatsPTM (2022)", "10.1016/j.mcpro.2022.100477"),
    ("MSstatsLiP / LiP-MS protocol (2023)", "10.1038/s41596-022-00771-x"),
    ("MSstatsShiny (2023)", "10.1021/acs.jproteome.2c00603"),
    ("MSstatsResponse (2026, preprint)", "10.64898/2026.03.09.710598"),
]

OPENALEX_URL = "https://api.openalex.org/works/https://doi.org/{doi}"


def fetch_cited_by_count(doi: str):
    req = urllib.request.Request(
        OPENALEX_URL.format(doi=doi),
        headers={"User-Agent": "msstats-readme-stats-bot (mailto:none)"},
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as response:
            data = json.load(response)
    except Exception:
        return None
    return data.get("cited_by_count")


def format_table(today: datetime.date) -> str:
    lines = [
        START_MARKER,
        "",
        "_Citation counts from [OpenAlex](https://openalex.org), updated "
        "monthly. Google Scholar counts are typically higher, since Scholar "
        "indexes a broader range of sources (theses, gray literature, etc.); "
        f"OpenAlex has a free, stable, official API. Last updated {today.isoformat()} (UTC)._",
        "",
        "| Paper | Citations |",
        "| --- | --- |",
    ]
    total = 0
    any_missing = False
    for label, doi in PAPERS:
        count = fetch_cited_by_count(doi)
        if count is None:
            any_missing = True
            cell = "n/a"
        else:
            total += count
            cell = f"{count:,}"
        lines.append(f"| [{label}](https://doi.org/{doi}) | {cell} |")

    total_cell = f"{total:,}" + ("+" if any_missing else "")
    lines.append(f"| **Total** | **{total_cell}** |")
    lines.append("")
    lines.append(END_MARKER)
    return "\n".join(lines)


def main() -> None:
    today = datetime.datetime.now(datetime.timezone.utc).date()
    new_block = format_table(today)

    text = README_PATH.read_text()
    if START_MARKER not in text or END_MARKER not in text:
        raise SystemExit(
            f"Could not find {START_MARKER} / {END_MARKER} markers in README.md"
        )
    before, rest = text.split(START_MARKER, 1)
    _old_block, after = rest.split(END_MARKER, 1)
    updated_text = before + new_block + after
    README_PATH.write_text(updated_text)


if __name__ == "__main__":
    main()
