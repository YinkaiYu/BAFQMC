#!/usr/bin/env python3
"""Build both documentation languages from the Markdown stored on GitHub."""
from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
REPO = "https://github.com/YinkaiYu/BAFQMC"
ASSETS = ROOT / "docs/assets"


def prepare(markdown, source, *, chinese=False):
    destination = {
        p.resolve(): (Path("index.md") if p == ROOT / "README.md" else
                      p.relative_to(ROOT / "docs/zh" if chinese else ROOT))
        for p in markdown
    }
    assets_target = Path("assets" if chinese else "docs/assets")
    if source.exists():
        shutil.rmtree(source)
    source.mkdir(parents=True)

    def rewrite_link(match, path, target):
        url = match.group(2)
        parsed = urlsplit(url)
        if parsed.scheme or parsed.netloc or not parsed.path or url.startswith("/"):
            return match.group(0)
        local = (path.parent / unquote(parsed.path)).resolve()
        fragment = "#" + parsed.fragment if parsed.fragment else ""
        if not chinese and local == ROOT / "README_zh-CN.md":
            replacement = "/BAFQMC/zh/" + fragment
        elif chinese and local == ROOT / "README.md":
            replacement = "/BAFQMC/" + fragment
        elif local in destination:
            replacement = Path(os.path.relpath(destination[local], target.parent)).as_posix() + fragment
        elif local.is_relative_to(ASSETS) and local.is_file():
            asset = assets_target / local.relative_to(ASSETS)
            replacement = Path(os.path.relpath(asset, target.parent)).as_posix() + fragment
        elif local.exists() and local.is_relative_to(ROOT):
            kind = "tree" if local.is_dir() else "blob"
            replacement = f"{REPO}/{kind}/main/{local.relative_to(ROOT).as_posix()}" + fragment
        else:
            raise ValueError(f"unresolved documentation link: {path.relative_to(ROOT)} -> {url}")
        return match.group(1) + replacement + match.group(3)

    for path in markdown:
        target = destination[path.resolve()]
        content = path.read_text().replace('<div align="center">', '<div align="center" markdown="1">')
        # GitHub math fences preserve TeX punctuation. Arithmatex also needs
        # blank lines around displays, even when a source fence touches prose.
        content = re.sub(r"(?m)^```math[ \t]*\n([\s\S]*?)^```[ \t]*$",
                         lambda m: "\n$$\n" + m.group(1) + "$$\n", content)
        content = re.sub(r"\$`([^`]+)`\$", lambda m: "$" + m.group(1) + "$", content)
        content = content.replace('<p align="center">', '<p align="center" markdown="1">')
        content = re.sub(r'<img src="([^"]+)" alt="([^"]+)" width="([^"]+)">',
                         r'![\2](\1){ width="\3" }', content)
        content = re.sub(r"(\]\()([^\s)]+)(\))",
                         lambda match: rewrite_link(match, path, target), content)
        out = source / target
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(content)
    shutil.copytree(ASSETS, source / assets_target)


def main():
    markdown = list(ROOT.glob("*.md"))
    for directory in ("docs", "benchmarks/paper", "src"):
        markdown.extend(p for p in (ROOT / directory).rglob("*.md")
                        if not {"build", "output", "data", "runs", "zh"}.intersection(p.relative_to(ROOT).parts))
    prepare(markdown, ROOT / ".build/docs-source")
    prepare(list((ROOT / "docs/zh").rglob("*.md")), ROOT / ".build/docs-zh-source", chinese=True)
    # English cleans the parent output first; Chinese is then built into /zh/.
    for config in ("mkdocs.yml", "mkdocs-zh.yml"):
        subprocess.run([sys.executable, "-m", "mkdocs", "build", "--strict",
                        "--config-file", str(ROOT / config)], cwd=ROOT, check=True)
    print(f"Documentation: {ROOT / '.build/site/index.html'} (English), "
          f"{ROOT / '.build/site/zh/index.html'} (简体中文)")


if __name__ == "__main__":
    main()
