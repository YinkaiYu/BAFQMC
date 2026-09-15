#!/usr/bin/env python3
"""Build the website from the same Markdown that readers and agents use on GitHub."""
from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / ".build/docs-source"
REPO = "https://github.com/YinkaiYu/BAFQMC"


def main():
    markdown = list(ROOT.glob("*.md"))
    for directory in ("docs", "benchmarks/paper", "src"):
        markdown.extend(p for p in (ROOT / directory).rglob("*.md") if not {"build", "output", "data", "runs"}.intersection(p.relative_to(ROOT).parts))
    destination = {p.resolve(): Path("index.md") if p.name == "README.md" and p.parent == ROOT else p.relative_to(ROOT) for p in markdown}
    if SOURCE.exists():
        shutil.rmtree(SOURCE)
    SOURCE.mkdir(parents=True)

    def rewrite_link(match, path, target):
        url = match.group(2)
        parsed = urlsplit(url)
        if parsed.scheme or parsed.netloc or not parsed.path or url.startswith("/"):
            return match.group(0)
        local = (path.parent / unquote(parsed.path)).resolve()
        fragment = "#" + parsed.fragment if parsed.fragment else ""
        if local in destination:
            replacement = Path(os.path.relpath(destination[local], target.parent)).as_posix() + fragment
        elif local.is_relative_to(ROOT / "docs/assets"):
            replacement = Path(os.path.relpath(local.relative_to(ROOT), target.parent)).as_posix() + fragment
        elif local.exists() and local.is_relative_to(ROOT):
            kind = "tree" if local.is_dir() else "blob"
            replacement = f"{REPO}/{kind}/main/{local.relative_to(ROOT).as_posix()}" + fragment
        else:
            raise ValueError(f"unresolved documentation link: {path.relative_to(ROOT)} -> {url}")
        return match.group(1) + replacement + match.group(3)

    for path in markdown:
        target = destination[path.resolve()]
        content = path.read_text().replace('<div align="center">', '<div align="center" markdown="1">')
        content = content.replace('<p align="center">', '<p align="center" markdown="1">')
        content = re.sub(r'<img src="([^"]+)" alt="([^"]+)" width="([^"]+)">', r'![\2](\1){ width="\3" }', content)
        content = re.sub(r"(\]\()([^\s)]+)(\))", lambda match: rewrite_link(match, path, target), content)
        out = SOURCE / target
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(content)
    shutil.copytree(ROOT / "docs/assets", SOURCE / "docs/assets")
    subprocess.run([sys.executable, "-m", "mkdocs", "build", "--strict", "--config-file", str(ROOT / "mkdocs.yml")], cwd=ROOT, check=True)
    print(f"Documentation: {ROOT / '.build/site/index.html'}")


if __name__ == "__main__":
    main()
