#!/usr/bin/env python3
"""Refresh the download-statistics table in README.md from Bioconductor's
per-package download stats, and write it between the
DOWNLOAD-STATS:START/END markers.

Uses only the standard library so it can run on a bare GitHub Actions
runner with no extra setup.
"""
import datetime
import pathlib
import urllib.request

README_PATH = pathlib.Path(__file__).resolve().parent.parent / "README.md"
START_MARKER = "<!-- DOWNLOAD-STATS:START -->"
END_MARKER = "<!-- DOWNLOAD-STATS:END -->"
MONTHS_TO_AVERAGE = 6

# (package name, display label) in the order they should appear in the table.
PACKAGES = [
    ("MSstats", "MSstats"),
    ("MSstatsTMT", "MSstatsTMT"),
    ("MSstatsPTM", "MSstatsPTM"),
    ("MSstatsLiP", "MSstatsLiP"),
    ("MSstatsBig", "MSstatsBig"),
    ("MSstatsShiny", "MSstatsShiny"),
    ("MSstatsBioNet", "MSstatsBioNet"),
    ("MSstatsResponse", "MSstatsResponse"),
    ("MSstatsConvert", "MSstatsConvert"),
]

MONTH_NUMBERS = {
    "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
    "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12,
}

STATS_URL = "https://bioconductor.org/packages/stats/bioc/{pkg}/{pkg}_stats.tab"


def fetch_stats(pkg: str) -> str:
    req = urllib.request.Request(
        STATS_URL.format(pkg=pkg),
        headers={"User-Agent": "msstats-readme-stats-bot/1.0"},
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return response.read().decode("utf-8")


def completed_months(raw_tab: str, today: datetime.date):
    """Yield (year, month_num, downloads) for months strictly before the
    current, still-in-progress month, most recent first."""
    rows = []
    for line in raw_tab.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) != 4:
            continue
        year_str, month_str, _ips, downloads_str = parts
        if month_str == "all":
            continue
        year = int(year_str)
        month = MONTH_NUMBERS[month_str]
        if (year, month) >= (today.year, today.month):
            continue
        rows.append((year, month, int(downloads_str)))
    rows.sort(reverse=True)
    return rows


def average_recent_downloads(pkg: str, today: datetime.date):
    raw_tab = fetch_stats(pkg)
    rows = completed_months(raw_tab, today)[:MONTHS_TO_AVERAGE]
    if not rows:
        return None, 0
    avg = round(sum(downloads for _, _, downloads in rows) / len(rows))
    return avg, len(rows)


def format_table(today: datetime.date) -> str:
    lines = [
        START_MARKER,
        "",
        f"_Average monthly downloads over the last {MONTHS_TO_AVERAGE} complete "
        f"months, computed directly from Bioconductor's download logs. Last "
        f"updated {today.isoformat()} (UTC)._",
        "",
        "| Package | Avg. monthly downloads |",
        "| --- | --- |",
    ]
    for pkg, label in PACKAGES:
        avg, n_months = average_recent_downloads(pkg, today)
        if avg is None:
            cell = "n/a (no completed months yet)"
        elif n_months < MONTHS_TO_AVERAGE:
            cell = f"{avg:,} (avg. over {n_months} available month{'s' if n_months != 1 else ''})"
        else:
            cell = f"{avg:,}"
        lines.append(
            f"| [{label}](https://bioconductor.org/packages/stats/bioc/{pkg}/) | {cell} |"
        )
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
