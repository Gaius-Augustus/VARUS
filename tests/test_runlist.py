"""Tests for varus.runlist.

We mock both Entrez.esearch and Entrez.esummary so the tests run offline. The
XML snippets below are minimal but structurally faithful to what the SRA
esummary endpoint returns.
"""

from __future__ import annotations

from io import BytesIO
from pathlib import Path

import pytest

from varus import runlist


# A two-DocSum esummary page: one paired Illumina run, one unpaired colorspace
# run on ABI_SOLID. The Run line uses the same attribute layout that SRA
# esummary actually returns inside the <Item Name="Runs"> CDATA.
ESUMMARY_XML = """<?xml version="1.0" ?>
<eSummaryResult>
<DocSum>
  <Id>1</Id>
  <Item Name="ExpXml" Type="String">
    &lt;LIBRARY_LAYOUT&gt;&lt;PAIRED/&gt;&lt;/LIBRARY_LAYOUT&gt;
    &lt;Instrument ILLUMINA="HiSeq 2500"/&gt;
  </Item>
  <Item Name="Runs" Type="String">
    &lt;Run acc="SRR000001" total_spots="1000" total_bases="100000"/&gt;
  </Item>
</DocSum>
<DocSum>
  <Id>2</Id>
  <Item Name="ExpXml" Type="String">
    &lt;LIBRARY_LAYOUT&gt;&lt;SINGLE/&gt;&lt;/LIBRARY_LAYOUT&gt;
    &lt;Instrument ABI_SOLID="AB SOLiD 4"/&gt;
  </Item>
  <Item Name="Runs" Type="String">
    &lt;Run acc="SRR000002" total_spots="500" total_bases="20000"/&gt;
  </Item>
</DocSum>
</eSummaryResult>
"""


def _fake_esearch_handle():
    return BytesIO(
        b"""<?xml version="1.0" ?>
<eSearchResult>
<Count>2</Count>
<RetMax>0</RetMax>
<RetStart>0</RetStart>
<QueryKey>1</QueryKey>
<WebEnv>FAKE_WEBENV</WebEnv>
</eSearchResult>"""
    )


def _fake_esummary_handle():
    return BytesIO(ESUMMARY_XML.encode("utf-8"))


def test_parse_xml_page_extracts_runs():
    records = list(runlist._parse_xml_page(ESUMMARY_XML))
    assert len(records) == 2
    paired = next(r for r in records if r.accession == "SRR000001")
    color = next(r for r in records if r.accession == "SRR000002")

    assert paired.total_spots == 1000
    assert paired.total_bases == 100000
    assert paired.avg_len == 100.0
    assert paired.paired is True
    assert paired.colorspace is False
    assert paired.platform == "ILLUMINA"

    assert color.paired is False
    assert color.colorspace is True
    assert color.avg_len == 40.0
    assert color.platform == "ABI_SOLID"


def test_parse_xml_page_detects_pacbio_platform():
    xml = ESUMMARY_XML.replace(
        'Instrument ILLUMINA="HiSeq 2500"',
        'Instrument PACBIO_SMRT="PacBio Sequel II"',
    )
    records = list(runlist._parse_xml_page(xml))
    pb = next(r for r in records if r.accession == "SRR000001")
    assert pb.platform == "PACBIO_SMRT"


def test_parse_xml_page_detects_ont_platform():
    xml = ESUMMARY_XML.replace(
        'Instrument ILLUMINA="HiSeq 2500"',
        'Instrument OXFORD_NANOPORE="PromethION"',
    )
    records = list(runlist._parse_xml_page(xml))
    ont = next(r for r in records if r.accession == "SRR000001")
    assert ont.platform == "OXFORD_NANOPORE"


def test_parse_xml_page_skips_blank_spots():
    xml = ESUMMARY_XML.replace('total_spots="1000"', 'total_spots=""')
    records = list(runlist._parse_xml_page(xml))
    assert {r.accession for r in records} == {"SRR000002"}


def test_fetch_runlist_writes_expected_tsv(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(
        runlist.Entrez, "esearch",
        lambda **kw: _fake_esearch_handle(),
    )
    monkeypatch.setattr(
        runlist.Entrez, "esummary",
        lambda **kw: _fake_esummary_handle(),
    )
    # Entrez.read on the esearch handle parses real XML, so let it run.

    out = runlist.fetch_runlist(
        species="Foo bar",
        outdir=tmp_path,
        email="test@example.org",
    )
    assert out == tmp_path / "Runlist.tsv"

    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("@Run_acc")
    assert "platform" in lines[0]
    body = [ln.split("\t") for ln in lines[1:]]
    accs = [row[0] for row in body]
    assert accs == ["SRR000001", "SRR000002"]
    # 7 columns now (added platform after color_space).
    assert all(len(row) == 7 for row in body)
    # Platform column populated from the SRA <Instrument> tag.
    platforms = [row[6] for row in body]
    assert platforms == ["ILLUMINA", "ABI_SOLID"]


def test_fetch_runlist_paired_only(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(
        runlist.Entrez, "esearch", lambda **kw: _fake_esearch_handle()
    )
    monkeypatch.setattr(
        runlist.Entrez, "esummary", lambda **kw: _fake_esummary_handle()
    )

    out = runlist.fetch_runlist(
        species="Foo bar",
        outdir=tmp_path,
        paired_only=True,
        email="test@example.org",
    )
    body = out.read_text().splitlines()[1:]
    accs = [ln.split("\t")[0] for ln in body]
    assert accs == ["SRR000001"]


def test_fetch_runlist_max_runs(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(
        runlist.Entrez, "esearch", lambda **kw: _fake_esearch_handle()
    )
    monkeypatch.setattr(
        runlist.Entrez, "esummary", lambda **kw: _fake_esummary_handle()
    )

    out = runlist.fetch_runlist(
        species="Foo bar",
        outdir=tmp_path,
        max_runs=1,
        email="test@example.org",
    )
    body = out.read_text().splitlines()[1:]
    assert len(body) == 1


def test_fetch_runlist_longreads_adds_platform_filter(tmp_path: Path, monkeypatch):
    """`--longreads` must restrict the Entrez term to PacBio + ONT platforms."""
    captured: dict = {}

    def _capturing_esearch(**kw):
        captured["term"] = kw.get("term")
        return _fake_esearch_handle()

    monkeypatch.setattr(runlist.Entrez, "esearch", _capturing_esearch)
    monkeypatch.setattr(
        runlist.Entrez, "esummary", lambda **kw: _fake_esummary_handle()
    )

    runlist.fetch_runlist(
        species="Foo bar",
        outdir=tmp_path,
        longreads=True,
        email="test@example.org",
    )
    assert "PACBIO_SMRT[Platform]" in captured["term"]
    assert "OXFORD_NANOPORE[Platform]" in captured["term"]


def test_fetch_runlist_default_omits_platform_filter(tmp_path: Path, monkeypatch):
    captured: dict = {}

    def _capturing_esearch(**kw):
        captured["term"] = kw.get("term")
        return _fake_esearch_handle()

    monkeypatch.setattr(runlist.Entrez, "esearch", _capturing_esearch)
    monkeypatch.setattr(
        runlist.Entrez, "esummary", lambda **kw: _fake_esummary_handle()
    )
    runlist.fetch_runlist(
        species="Foo bar", outdir=tmp_path, email="test@example.org",
    )
    assert "Platform" not in captured["term"]


def test_fetch_runlist_raises_when_no_runs(tmp_path: Path, monkeypatch):
    def _empty_esearch(**kw):
        return BytesIO(
            b"""<?xml version="1.0" ?>
<eSearchResult><Count>0</Count><WebEnv>X</WebEnv><QueryKey>1</QueryKey></eSearchResult>"""
        )

    monkeypatch.setattr(runlist.Entrez, "esearch", _empty_esearch)
    with pytest.raises(RuntimeError, match="No SRA RNA-seq runs found"):
        runlist.fetch_runlist(
            species="Nope nope",
            outdir=tmp_path,
            email="test@example.org",
        )
