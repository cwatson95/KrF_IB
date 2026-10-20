"""Integrity checks for the repository's browsable codebase guide."""

from __future__ import annotations

import unittest
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
GUIDE_PATH = REPOSITORY_ROOT / "docs" / "index.html"


class _GuideParser(HTMLParser):
    """Collect document structure, identifiers, and navigation targets."""

    def __init__(self) -> None:
        super().__init__()
        self.ids: list[str] = []
        self.hrefs: list[str] = []
        self.elements: Counter[str] = Counter()

    def handle_starttag(
        self,
        tag: str,
        attrs: list[tuple[str, str | None]],
    ) -> None:
        self.elements[tag] += 1
        attributes = dict(attrs)
        if identifier := attributes.get("id"):
            self.ids.append(identifier)
        if href := attributes.get("href"):
            self.hrefs.append(href)


class TestCodebaseGuide(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.parser = _GuideParser()
        cls.parser.feed(GUIDE_PATH.read_text(encoding="utf-8"))
        cls.parser.close()

    def test_document_has_one_complete_html_structure(self) -> None:
        for element in ("html", "head", "title", "body", "main"):
            with self.subTest(element=element):
                self.assertEqual(self.parser.elements[element], 1)

    def test_ids_are_unique(self) -> None:
        duplicates = [
            identifier
            for identifier, count in Counter(self.parser.ids).items()
            if count > 1
        ]
        self.assertEqual(duplicates, [])

    def test_fragment_links_resolve(self) -> None:
        identifiers = set(self.parser.ids)
        missing = sorted(
            href for href in self.parser.hrefs if href.startswith("#") and href[1:] not in identifiers
        )
        self.assertEqual(missing, [])

    def test_local_links_resolve_inside_repository(self) -> None:
        unresolved: list[str] = []
        escaping: list[str] = []
        for href in self.parser.hrefs:
            parsed = urlsplit(href)
            if parsed.scheme or not parsed.path:
                continue
            target = (GUIDE_PATH.parent / unquote(parsed.path)).resolve()
            if not target.is_relative_to(REPOSITORY_ROOT):
                escaping.append(href)
            elif not target.exists():
                unresolved.append(href)
        self.assertEqual(escaping, [])
        self.assertEqual(unresolved, [])

    def test_readme_links_to_full_guide(self) -> None:
        readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn("](docs/index.html)", readme)


if __name__ == "__main__":
    unittest.main()
