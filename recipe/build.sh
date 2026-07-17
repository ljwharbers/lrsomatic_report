#!/bin/bash
set -euo pipefail

target="$PREFIX/share/lrsomatic_report"
mkdir -p "$target"
cp -r bin R templates assets "$target/"
chmod +x "$target/bin/render_report.R"

mkdir -p "$PREFIX/bin"
ln -s ../share/lrsomatic_report/bin/render_report.R "$PREFIX/bin/render_report.R"
