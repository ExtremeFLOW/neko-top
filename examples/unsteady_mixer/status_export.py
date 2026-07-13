#!/usr/bin/env python3
"""Export a status notebook to markdown with embedded image files.

This script mirrors the notebook export cell behavior in a pure Python entrypoint:
1. load notebook,
2. drop export cells to avoid recursion,
3. execute notebook,
4. export markdown with code inputs hidden,
5. store output images in a dedicated folder and rewrite markdown links.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=
        "Execute a Jupyter notebook and export markdown with images.", )
    parser.add_argument(
        "--notebook",
        default="status.ipynb",
        help="Notebook file to export (default: status.ipynb)",
    )
    parser.add_argument(
        "--markdown",
        default="status.md",
        help="Output markdown filename (default: status.md)",
    )
    parser.add_argument(
        "--image-dir",
        default="status",
        help="Directory (relative to notebook directory) for extracted images.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Per-cell execution timeout in seconds (default: 600)",
    )
    return parser.parse_args()


def _replace_markdown_image_links(markdown: str, filename: str,
                                  image_dir: str) -> str:
    """Rewrite markdown links for extracted files to point under image_dir."""
    escaped = re.escape(filename)
    return re.sub(rf"\(({escaped})\)", f"({image_dir}/{filename})", markdown)


def main() -> int:
    args = _parse_args()

    script_dir = Path(__file__).resolve().parent
    notebook_arg = Path(args.notebook)
    notebook_path = (notebook_arg.resolve() if notebook_arg.is_absolute() else
                     (script_dir / notebook_arg).resolve())
    if not notebook_path.exists():
        print(f"Notebook not found: {notebook_path}", file=sys.stderr)
        return 1

    workdir = notebook_path.parent
    markdown_path = workdir / args.markdown
    image_dir = workdir / args.image_dir
    image_dir.mkdir(parents=True, exist_ok=True)

    try:
        import nbformat
        from nbconvert import MarkdownExporter
        from nbconvert.preprocessors import ExecutePreprocessor
    except ImportError as exc:
        print(
            "Missing dependency for notebook export. "
            "Install with: pip install nbformat nbconvert",
            file=sys.stderr,
        )
        print(f"Import error: {exc}", file=sys.stderr)
        return 2

    with notebook_path.open(encoding="utf-8") as f:
        notebook_node = nbformat.read(f, as_version=4)

    notebook_node.cells = [
        cell for cell in notebook_node.cells if not str(cell.get(
            "source", "")).lstrip().startswith("# Export Notebook as")
    ]
    notebook_node.metadata.title = "Status of optimization"

    executor = ExecutePreprocessor(timeout=args.timeout, kernel_name="python3")
    executor.preprocess(notebook_node, {"metadata": {"path": str(workdir)}})

    exporter = MarkdownExporter(exclude_input=True, exclude_input_prompt=True)
    markdown, resources = exporter.from_notebook_node(notebook_node)

    for filename, data in resources.get("outputs", {}).items():
        out_path = image_dir / filename
        with out_path.open("wb") as f:
            f.write(data)
        markdown = _replace_markdown_image_links(markdown, filename,
                                                 args.image_dir)

    with markdown_path.open("w", encoding="utf-8") as f:
        f.write(markdown)

    print(f"Exported markdown: {markdown_path}")
    print(f"Exported images dir: {image_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
