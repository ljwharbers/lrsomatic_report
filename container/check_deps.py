#!/usr/bin/env python3
"""Fail if container/environment.yml and recipe/meta.yaml name different packages.

The image is built from environment.yml and the bioconda package from meta.yaml, so a
dependency added to one and not the other produces two builds of "the same" tool that do
not carry the same software. Versions are deliberately not compared: the recipe keeps loose
pins because it is a noarch package, while the image pins what it was tested against.
"""

import pathlib
import re
import sys

DEP = re.compile(r"\s*-\s*(?:[a-z0-9-]+::)?([A-Za-z0-9._-]+)")


def container_deps(path):
    lines = pathlib.Path(path).read_text().splitlines()
    start = lines.index("dependencies:")
    names = set()
    for line in lines[start + 1 :]:
        if line.strip() and not line.startswith((" ", "\t", "-")):
            break
        m = DEP.match(line)
        if m:
            names.add(m.group(1))
    return names


def recipe_run_deps(path):
    lines = pathlib.Path(path).read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip() == "run:")
    indent = len(lines[start]) - len(lines[start].lstrip())
    names = set()
    for line in lines[start + 1 :]:
        if not line.strip():
            continue
        # Leaving the run: block (back to the same or lower indent, e.g. the next section)
        if len(line) - len(line.lstrip()) <= indent and not line.lstrip().startswith("-"):
            break
        m = DEP.match(line)
        if m:
            names.add(m.group(1))
    return names


def main():
    root = pathlib.Path(__file__).resolve().parent.parent
    container = container_deps(root / "container/environment.yml")
    recipe = recipe_run_deps(root / "recipe/meta.yaml")
    if container != recipe:
        print("::error::container/environment.yml and recipe/meta.yaml disagree")
        print("  only in container/environment.yml:", sorted(container - recipe))
        print("  only in recipe/meta.yaml:         ", sorted(recipe - container))
        return 1
    print(f"dependency lists agree ({len(container)} packages)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
