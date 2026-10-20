"""LIGHT_Kr v0_4's two documents (LIGHT_Kr/docs): the accuracy overview
and the system map. Held here:

* every number the .tex files print is defined in the generated numbers
  files, and the numbers were generated from the model's current code
  (the code hash each file records);
* every citation has a bibliography entry, every entry is cited and
  carries a DOI, and the text names no local file or path: the documents
  stand alone (the owner's rule: only published sources with a DOI);
* with LIGHT_KINETIC_SLOW=1, both generators reproduce their files
  exactly (minutes from a current run cache, about half an hour without).

Skips when LIGHT_Kr/docs is not beside v0_4."""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parents[2] / "docs"
MAP_TEX = DOCS / "system_map_v0_4.tex"
OV_TEX = DOCS / "accuracy_overview_v0_4.tex"
MAP_NUM = DOCS / "system_map_v0_4_numbers.tex"
OV_NUM = DOCS / "overview_v0_4_numbers.tex"
needs_docs = pytest.mark.skipif(not OV_TEX.is_file(),
                                reason="LIGHT_Kr/docs is not beside v0_4")
SLOW = pytest.mark.skipif(not os.environ.get("LIGHT_KINETIC_SLOW"),
                          reason="minutes: set LIGHT_KINETIC_SLOW=1")


def _defined(path, macro):
    return set(re.findall(r"\\" + macro + r"def\{([^}]*)\}",
                          path.read_text()))


def _used(text, macro):
    """The keys a document prints, not a macro's own parameters (#1, \\k
    in a loop)."""
    return {k for k in re.findall(r"\\" + macro + r"\{([^}]*)\}", text)
            if not k.startswith(("#", "\\"))}


def _body(path):
    """The document's body without its comment lines and comments."""
    text = path.read_text()
    body = text.split(r"\begin{document}", 1)[1]
    lines = []
    for line in body.splitlines():
        line = re.sub(r"(?<!\\)%.*", "", line)
        if line.strip():
            lines.append(line)
    return "\n".join(lines)


@needs_docs
def test_every_number_printed_is_generated():
    wm = _defined(MAP_NUM, "wm")
    wg = _defined(OV_NUM, "wg")
    for tex in (MAP_TEX, OV_TEX):
        text = tex.read_text()
        missing = _used(text, "wm") - wm
        assert not missing, (tex.name, sorted(missing))
    missing = _used(OV_TEX.read_text(), "wg") - wg
    assert not missing, sorted(missing)
    for name in ("wgUncTable", "wgScan", "wgSwTable", "wgFacTable",
                 "wgVarTable", "wgMcTable", "wgSwWtr", "wgSwWbi", "wgSwAtr",
                 "wgSwAbi", "wgSwAHM", "wgSwDtr", "wgSwDbi", "wgSwWTdata",
                 "wgSwWHdata", "wgSwATdata", "wgSwAHdata", "wgSwDHdata",
                 # v0_4: what each change does, the closure against the
                 # kinetic EEDF, the electroluminescence score
                 "wgUndo", "wgTail", "wgEl"):
        assert "\\def\\" + name + "{" in OV_NUM.read_text()


@needs_docs
def test_numbers_come_from_the_current_code():
    sys.path.insert(0, str(DOCS))
    try:
        import kr_docs_v0_4 as kr_docs
    finally:
        sys.path.remove(str(DOCS))
    # v0_3: the stamp covers every module and data file, the swarm modules
    # the generators read live included (kr_docs_v0_4.full_md5)
    md5 = kr_docs.full_md5()[:12]
    for path, macro in ((MAP_NUM, "wm"), (OV_NUM, "wg")):
        stamp = re.search(r"\\" + macro + r"def\{md5\}\{([0-9a-f]+)\}",
                          path.read_text())
        assert stamp and stamp.group(1) == md5, (
            f"{path.name} was made from other code: rerun its generator")


@needs_docs
def test_citations_are_published_sources_with_a_doi():
    text = OV_TEX.read_text()
    cited = set()
    for keys in re.findall(r"\\cite\{([^}]*)\}", text):
        cited.update(k.strip() for k in keys.split(","))
    items = dict(re.findall(r"\\bibitem\{([^}]*)\}(.*)", text))
    assert cited == set(items), (sorted(cited - set(items)),
                                 sorted(set(items) - cited))
    for key, entry in items.items():
        assert re.search(r"doi:10\.\d{4,9}/\S+", entry), key


@needs_docs
@pytest.mark.parametrize("tex", [MAP_TEX, OV_TEX], ids=lambda p: p.name)
def test_the_text_names_no_local_file(tex):
    body = _body(tex)
    for pattern in (r"\w\.py\b", r"\.md\b", r"\.txt\b", r"\.csv\b",
                    r"\.pdf\b", r"\w/docs\b", r"docs/", r"v0\\?_\d/",
                    r"SOURCES", r"README", r"CHANGELOG", r"Kinetic_Modeling",
                    r"/Users/", r"~/", r"this repository"):
        hit = re.search(pattern, body)
        assert hit is None, (tex.name, pattern, body[max(0, hit.start()
                                                         - 60):hit.end()
                                                    + 20])


@needs_docs
@SLOW
@pytest.mark.parametrize("script", ["make_system_map_v0_4.py",
                                    "make_overview_v0_4.py"])
def test_generated_files_are_current(script):
    r = subprocess.run([sys.executable, str(DOCS / script), "--check"],
                       capture_output=True, text=True, cwd=DOCS)
    assert r.returncode == 0, r.stdout + r.stderr
