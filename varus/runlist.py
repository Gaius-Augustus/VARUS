"""Fetch the SRA RNA-seq run list for a species via NCBI Entrez.

Replaces ``legacy/RunListRetriever/RunListRetriever.pl``.

The Perl version did ``wget`` + regex on the raw XML. We use ``Bio.Entrez`` with
``usehistory=y`` and page through ``esummary`` 10 000 records at a time, matching
the page size the original tool used.

Output ``Runlist.tsv`` has the same columns as the legacy ``Runlist.txt`` so
existing downstream tooling keeps working::

    @Run_acc  total_spots  total_bases  avg_len  bool:paired  color_space
"""

from __future__ import annotations

import logging
import os
import random
import re
import time
import urllib.error
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator
from xml.etree import ElementTree as ET

from Bio import Entrez

log = logging.getLogger(__name__)

PAGE_SIZE = 10_000  # NCBI esummary cap per call


@dataclass(frozen=True)
class RunRecord:
    accession: str
    total_spots: int
    total_bases: int
    avg_len: float
    paired: bool
    colorspace: bool
    # SRA platform identifier (uppercase: e.g. ILLUMINA, PACBIO_SMRT,
    # OXFORD_NANOPORE, ABI_SOLID, ION_TORRENT, BGISEQ). Empty string when the
    # esummary XML lacked an <Instrument> tag we could parse.
    platform: str = ""


def _configure_entrez(email: str | None, api_key: str | None) -> None:
    Entrez.email = email or os.environ.get("NCBI_EMAIL", "varus@example.org")
    api_key = api_key or os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key


# Entrez SRA Platform[Platform] filter for long-read RNA-seq submissions.
# Covers the two long-read platforms in current production use; can be extended
# if a third major one (e.g. Element Bio) becomes common in SRA.
LONGREAD_PLATFORM_TERM = "(PACBIO_SMRT[Platform] OR OXFORD_NANOPORE[Platform])"


def _esearch_history(species: str, longreads: bool = False,
                     retries: int = 6) -> tuple[int, str, str]:
    """Run esearch with ``usehistory=y``; return (count, WebEnv, query_key).

    Parses the response with stdlib ElementTree rather than ``Entrez.read``
    because the latter requires a DTD reference in the XML.

    Retries on HTTP 429 (NCBI rate-limit) with jittered exponential backoff,
    so parallel pipeline invocations don't permanently drop species when they
    momentarily exceed the unauthenticated 3 req/s cap. A genuine empty
    result (count == 0) is NOT retried — it's a terminal "no data" signal.
    """
    term = f'"{species}"[orgn] AND biomol_rna[Prop]'
    if longreads:
        term += f" AND {LONGREAD_PLATFORM_TERM}"
    log.info("Entrez esearch term=%s", term)
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="sra", term=term, usehistory="y", retmax=0)
            try:
                data = handle.read()
            finally:
                handle.close()
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                wait = (2 ** attempt) + random.random()
                log.warning("Entrez esearch HTTP %s on attempt %d for %r; "
                            "retrying in %.1fs", e.code, attempt + 1, species, wait)
                time.sleep(wait)
                last_err = e
                continue
            raise
        if isinstance(data, bytes):
            data = data.decode("utf-8")
        root = ET.fromstring(data)
        count = int(root.findtext("Count", "0"))
        if count == 0:
            raise RuntimeError(f"No SRA RNA-seq runs found for species {species!r}")
        webenv = root.findtext("WebEnv") or ""
        query_key = root.findtext("QueryKey") or "1"
        return count, webenv, query_key
    raise RuntimeError(
        f"Entrez esearch failed after {retries} retries for {species!r}: {last_err}"
    )


def _esummary_page(webenv: str, query_key: str, retstart: int,
                   retmax: int = PAGE_SIZE, retries: int = 3) -> str:
    """Fetch one page of esummary XML as a string."""
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            handle = Entrez.esummary(
                db="sra",
                WebEnv=webenv,
                query_key=query_key,
                retstart=retstart,
                retmax=retmax,
            )
            try:
                data = handle.read()
            finally:
                handle.close()
            if isinstance(data, bytes):
                data = data.decode("utf-8")
            return data
        except Exception as e:  # network glitches, transient 5xx
            last_err = e
            wait = 2 ** attempt
            log.warning("esummary retstart=%d attempt %d failed: %s; "
                        "retrying in %ds", retstart, attempt + 1, e, wait)
            time.sleep(wait)
    raise RuntimeError(f"esummary failed after {retries} retries: {last_err}")


# The SRA esummary returns one DocSum per run. The Run accession, total_spots
# and total_bases live in an XML attribute string inside <Item Name="Runs">,
# not as proper child elements -- which is why the legacy code regexed it. We
# do the same thing here, but on a more forgiving regex than the Perl one.
_RUN_RE = re.compile(
    r'Run\s+acc="([^"]+)"\s+total_spots="(\d*)"\s+total_bases="(\d*)"',
)
_LAYOUT_RE = re.compile(r"LAYOUT.{0,32}PAIRED", re.DOTALL)
_COLORSPACE_RE = re.compile(r"Instrument\s+ABI_SOLID")
# Captures the SRA platform attribute on the <Instrument> tag, e.g.
# `<Instrument PACBIO_SMRT="PacBio Sequel II"/>` -> "PACBIO_SMRT". Encoded as
# `&lt;Instrument PACBIO_SMRT=...` in the esummary CDATA, so match either form.
_PLATFORM_RE = re.compile(r"Instrument\s+([A-Z][A-Z0-9_]+)\s*=")


def _parse_xml_page(xml: str) -> Iterator[RunRecord]:
    """Yield RunRecord per DocSum block found in this page.

    The DocSum XML for SRA contains ``<ExpXml>`` and ``<Runs>`` items whose
    *contents* are escaped XML. Rather than nested-parse, we split the page on
    DocSum boundaries and regex the inner block, which is what the legacy tool
    effectively did via line-streaming.
    """
    blocks = xml.split("<DocSum>")
    for block in blocks[1:]:
        end = block.find("</DocSum>")
        if end >= 0:
            block = block[:end]
        paired = bool(_LAYOUT_RE.search(block))
        colorspace = bool(_COLORSPACE_RE.search(block))
        platform_m = _PLATFORM_RE.search(block)
        platform = platform_m.group(1) if platform_m else ""
        for m in _RUN_RE.finditer(block):
            acc, spots_s, bases_s = m.group(1), m.group(2), m.group(3)
            if not spots_s or not bases_s:
                continue
            spots = int(spots_s)
            bases = int(bases_s)
            if spots <= 0:
                continue
            avg_len = round(100 * bases / spots) / 100.0
            yield RunRecord(
                accession=acc,
                total_spots=spots,
                total_bases=bases,
                avg_len=avg_len,
                paired=paired,
                colorspace=colorspace,
                platform=platform,
            )


def _iter_runs(species: str, longreads: bool = False) -> Iterator[RunRecord]:
    count, webenv, query_key = _esearch_history(species, longreads=longreads)
    log.info("Server has %d data sets for %s", count, species)
    retstart = 0
    while retstart < count:
        log.info("esummary page retstart=%d", retstart)
        xml = _esummary_page(webenv, query_key, retstart)
        yield from _parse_xml_page(xml)
        retstart += PAGE_SIZE


def write_runlist(records: Iterable[RunRecord], path: Path) -> int:
    n = 0
    with path.open("w", encoding="utf-8") as f:
        f.write(
            "@Run_acc\ttotal_spots\ttotal_bases\tavg_len\tbool:paired\t"
            "color_space\tplatform\n"
        )
        for r in records:
            f.write(
                f"{r.accession}\t{r.total_spots}\t{r.total_bases}\t"
                f"{r.avg_len}\t{int(r.paired)}\t{int(r.colorspace)}\t"
                f"{r.platform}\n"
            )
            n += 1
    return n


def fetch_runlist(
    species: str,
    outdir: Path,
    max_runs: int = 0,
    paired_only: bool = False,
    longreads: bool = False,
    email: str | None = None,
    api_key: str | None = None,
) -> Path:
    """Query SRA, filter, and write ``<outdir>/Runlist.tsv``.

    When ``longreads`` is True the Entrez term is restricted to PacBio SMRT and
    Oxford Nanopore platforms; without this filter SRA returns mostly Illumina
    short-read runs because they dominate the archive.

    Returns the path to the written file.
    """
    _configure_entrez(email, api_key)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "Runlist.tsv"

    def _filter(it: Iterator[RunRecord]) -> Iterator[RunRecord]:
        kept = 0
        for r in it:
            if paired_only and not r.paired:
                continue
            yield r
            kept += 1
            if max_runs and kept >= max_runs:
                break

    n = write_runlist(_filter(_iter_runs(species, longreads=longreads)), out)
    log.info("Wrote %d runs to %s", n, out)
    return out
