#!/usr/bin/env bash
# Verify recipe/meta.yaml's sha256 against the release tarball it names.
#
# The recipe has shipped a wrong digest more than once: the file at v1.3.2 still said 1.3.1,
# and main once paired version 1.4.0 with the v1.3.1 tarball's hash. Nothing caught either,
# because the digest can only be recorded after the tag is cut -- so the release commit and
# the commit that fixes it are always separate, and the window between them is invisible.
#
# This is a no-op until that tag exists, which is what makes it safe to run on every PR:
# the release PR skips, and the PR that records the digest is checked.
set -euo pipefail

cd "$(dirname "$0")/.."

version=$(sed -n 's/^{%[[:space:]]*set version[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' recipe/meta.yaml)
recorded=$(sed -n 's/^[[:space:]]*sha256:[[:space:]]*\([0-9a-f]\{64\}\).*/\1/p' recipe/meta.yaml)

[ -n "$version" ]  || { echo "could not read the version from recipe/meta.yaml"; exit 1; }
[ -n "$recorded" ] || { echo "could not read the sha256 from recipe/meta.yaml"; exit 1; }

if [ "$version" != "$(cat VERSION)" ]; then
    echo "::error::recipe/meta.yaml says $version but VERSION says $(cat VERSION)"
    exit 1
fi

url="https://github.com/ljwharbers/lrsomatic_report/archive/refs/tags/v${version}.tar.gz"
if ! curl -fsSLI "$url" > /dev/null 2>&1; then
    echo "v${version} is not tagged yet, so there is no tarball to hash -- skipping"
    exit 0
fi

actual=$(curl -fsSL "$url" | sha256sum | cut -d' ' -f1)
if [ "$recorded" != "$actual" ]; then
    echo "::error::recipe/meta.yaml sha256 does not match the v${version} tarball"
    echo "  recorded: $recorded"
    echo "  actual:   $actual"
    exit 1
fi

echo "recipe/meta.yaml sha256 matches the v${version} tarball"
