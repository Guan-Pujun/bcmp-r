"""Check local documentation links and source-archive contents without deployment."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from html.parser import HTMLParser
from pathlib import Path
import re
import tarfile
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids = set()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.add(attrs["id"])
        for key in ("href", "src", "data-theme-src-light", "data-theme-src-dark"):
            if attrs.get(key):
                self.links.append(attrs[key])
        if attrs.get("srcset"):
            self.links.extend(item.strip().split()[0] for item in attrs["srcset"].split(","))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("docs", "tarball"))
    parser.add_argument("--tarball", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    errors, external = [], set()
    report = {"mode": args.mode, "errors": errors}
    if args.mode == "docs":
        pages = {p.resolve(): Page(p.read_text()) for p in (root / "docs").rglob("*.html")}
        if not pages:
            errors.append("No manual pages found")
        links = [(page, link) for page, data in pages.items() for link in data.links]
        readme = root / "README.md"
        links.extend((readme, link) for link in Page(readme.read_text()).links)
        links.extend((readme, link) for link in re.findall(r"\]\(([^)]+)\)", readme.read_text()))
        for css in (root / "docs").rglob("*.css"):
            links.extend((css, link) for link in re.findall(r"""url\(["']?([^)'" ]+)""", css.read_text()))
        for source, link in links:
            url = urlsplit(link)
            if url.scheme or url.netloc:
                if url.hostname in ("localhost", "127.0.0.1"):
                    errors.append(f"Local-only URL: {source.relative_to(root)} -> {link}")
                elif url.scheme != "data":
                    external.add(link)
                continue
            target = (source.parent / unquote(url.path)).resolve() if url.path else source
            if not target.is_relative_to(root) or not target.is_file():
                errors.append(f"Missing/unsafe local target: {source.relative_to(root)} -> {link}")
            elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
                errors.append(f"Missing anchor: {source.relative_to(root)} -> {link}")
        report.update(pages=len(pages), external_links=sorted(external), external_links_checked=False)
    else:
        if args.tarball is None:
            parser.error("--tarball is required")
        required = {"DESCRIPTION", "NAMESPACE"}
        for directory in ("R", "src", "man", "inst", "tests"):
            required.update(p.relative_to(root).as_posix() for p in (root / directory).rglob("*")
                            if p.is_file() and not p.name.endswith((".o", ".so", ".dll", ".orig")))
        files = set()
        with tarfile.open(args.tarball) as archive:
            for member in archive.getmembers():
                parts = Path(member.name).parts
                if not parts or parts[0] != "bcmp" or ".." in parts or member.issym() or member.islnk():
                    errors.append(f"Unsafe archive member: {member.name}")
                    continue
                if not member.isfile():
                    continue
                relative = "/".join(parts[1:])
                files.add(relative)
                if relative not in required:
                    errors.append(f"Unexpected source-package file: {relative}")
                    continue
                expected = root / relative
                if relative != "DESCRIPTION" and archive.extractfile(member).read() != expected.read_bytes():
                    errors.append(f"Archive/source mismatch: {relative}")
        if required - files:
            errors.append(f"Missing source-package files: {sorted(required - files)}")
        report.update(files=sorted(files), sha256=hashlib.sha256(args.tarball.read_bytes()).hexdigest())
        golden = root / "tests/testthat/fixtures/r-v1/diabetic-kidney-lite"
        with (golden / "checksums.csv").open() as stream:
            for row in csv.DictReader(stream):
                data = (golden / row["file"]).read_bytes()
                if hashlib.sha256(data).hexdigest() != row["sha256"]:
                    errors.append(f"Frozen checksum mismatch: {row['file']}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
