#!/usr/bin/env bash
# Deferred-start helper: waits until a given wall-clock time, then runs a
# long task in the repo. Fill in the placeholders before use.
#
# Kept as a template only -- do not commit real sample IDs, report filenames,
# or absolute paths under a user's home directory. This repository is public.

REPO=/path/to/lrsomatic_report
START=17:30

sleep $(( $(date -d "$START today" +%s) - $(date +%s) )) \
  && cd "$REPO" \
  && claude -p 'Execute the plan at <plan-path>. Render SAMPLE_ID at the end
      (per the Verification section) and save the HTML outside the repo.' \
      > ~/claude-run-$(date +%Y%m%d-%H%M).log 2>&1
