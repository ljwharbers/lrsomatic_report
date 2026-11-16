// Clinical overview: the key-findings table follows the gene-list chips.
//
// Inlined at the end of templates/per_sample.qmd by the key-findings-js chunk (same route as
// facet_filter.js), i.e. after every section's markup and after the panel-js-data globals.
// The markup is templates/sections/_overview.qmd: one <tr class="kf-row"> per finding with
//   data-gene   upper-case symbol      data-panels  comma-joined panel keys ("" = no panel)
//   data-tier   1 | 2 | 3              data-type    small variant | SV | CNA
// Rows are computed in R for EVERY embedded panel (R/key_findings.R); this script only shows
// and hides them. A row is shown when its panels intersect the ticked ones, or its gene is in
// a typed Custom… list. Nothing ticked ("all") shows no findings at all, but a notice with
// per-list totals that tick the chip when clicked: a findings table without a gene list would
// be every HIGH-impact variant in the genome, which is not a driver summary.
//
// window.registerPanelListener() calls the listener at least twice per change, so apply() is
// a pure function of (state, showLow) and safe to repeat.
//
// The predicates are exported for node (tests/js/test-key-findings.js); the DOM half only runs
// when a document exists.
(function () {
  "use strict";

  function asList(x) {
    if (!x) return [];
    if (Array.isArray(x)) return x;
    if (typeof x === "string") return x.split(",");
    return Array.from(x);
  }

  // Does one finding row match the panel state? row = {gene, panels: [...], tier}
  function kfRowMatches(row, state) {
    if (!state || state.all) return false;
    const panels = asList(state.panels);
    const rowPanels = asList(row.panels).filter(function (p) { return p.length > 0; });
    for (let i = 0; i < rowPanels.length; i++) {
      if (panels.indexOf(rowPanels[i]) >= 0) return true;
    }
    const custom = asList(state.customGenes).map(function (g) { return String(g).toUpperCase(); });
    return custom.indexOf(String(row.gene || "").toUpperCase()) >= 0;
  }

  // Visibility of every row plus the counts the widget prints. Tier-3 rows that match stay
  // hidden unless showLow; they are counted in nLow so the toggle can say how many.
  function kfSummarise(rows, state, showLow) {
    const all = !state || !!state.all;
    const visible = [];
    let matched = 0, nLow = 0, shown = 0;
    rows.forEach(function (row) {
      const m = kfRowMatches(row, state);
      const low = Number(row.tier) >= 3;
      if (m) matched++;
      if (m && low) nLow++;
      const v = m && (!low || !!showLow);
      if (v) shown++;
      visible.push(v);
    });
    return { visible: visible, matched: matched, nLow: nLow, shown: shown, all: all,
             empty: !all && matched === 0 };
  }

  // Human label of the active lists, for the print line
  function kfListLabel(state, labels) {
    if (!state || state.all) return "none (no gene list ticked)";
    const out = asList(state.panels).map(function (p) { return (labels && labels[p]) || p; });
    const custom = asList(state.customGenes);
    if (custom.length) out.push("Custom (" + custom.length + " gene" + (custom.length === 1 ? "" : "s") + ")");
    return out.length ? out.join(", ") : "none";
  }

  const api = { kfRowMatches: kfRowMatches, kfSummarise: kfSummarise, kfListLabel: kfListLabel };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  if (typeof window !== "undefined") window.keyFindings = api;
  if (typeof document === "undefined") return;

  const root = document.getElementById("kf-root");
  if (!root) return;
  const trs = Array.prototype.slice.call(root.querySelectorAll("tr.kf-row"));
  const rows = trs.map(function (tr) {
    return { gene: tr.getAttribute("data-gene") || "",
             panels: (tr.getAttribute("data-panels") || "").split(","),
             tier: Number(tr.getAttribute("data-tier")) };
  });
  const tableWrap = root.querySelector(".kf-table-wrap");
  const needList = root.querySelector(".kf-need-list");
  const empty = root.querySelector(".kf-empty");
  const toggle = root.querySelector(".kf-toggle-low");
  const count = root.querySelector(".kf-count");
  const printLine = root.querySelector(".kf-print-panels");
  const customNote = root.querySelector(".kf-custom-note");

  // Chip labels as the chip row prints them ("lymphoid" -> "Lymphoid")
  const labels = {};
  Array.prototype.forEach.call(document.querySelectorAll("input.panel-checkbox"), function (cb) {
    const lab = cb.closest ? cb.closest("label") : cb.parentNode;
    labels[cb.value] = lab ? lab.textContent.trim() : cb.value;
  });

  let state = { panels: [], customGenes: [], all: true };
  let showLow = false;

  function apply() {
    const s = kfSummarise(rows, state, showLow);
    trs.forEach(function (tr, i) { tr.hidden = !s.visible[i]; });
    if (tableWrap) tableWrap.hidden = s.all || s.matched === 0;
    if (needList) needList.hidden = !s.all;
    if (empty) empty.hidden = !s.empty;
    if (customNote) customNote.hidden = asList(state.customGenes).length === 0;
    if (toggle) {
      toggle.hidden = s.all || s.nLow === 0;
      toggle.setAttribute("aria-expanded", showLow ? "true" : "false");
      toggle.textContent = showLow
        ? "Hide " + s.nLow + " lower-confidence finding" + (s.nLow === 1 ? "" : "s")
        : "Show " + s.nLow + " lower-confidence finding" + (s.nLow === 1 ? "" : "s") + " (tier 3)";
    }
    if (count) {
      count.textContent = s.all ? "" :
        (s.shown + " of " + s.matched + " finding" + (s.matched === 1 ? "" : "s") + " shown");
    }
    if (printLine) printLine.textContent = "Gene list(s): " + kfListLabel(state, labels);
  }

  if (toggle) toggle.addEventListener("click", function () { showLow = !showLow; apply(); });

  // "Lymphoid (4)" hints in the empty state tick the chip itself, through its own change
  // handler (click() toggles the box and fires change), so tables and filters follow too
  Array.prototype.forEach.call(root.querySelectorAll(".kf-hint[data-panel]"), function (btn) {
    btn.addEventListener("click", function () {
      const key = btn.getAttribute("data-panel");
      Array.prototype.forEach.call(document.querySelectorAll("input.panel-checkbox"), function (cb) {
        if (cb.value === key && !cb.checked) cb.click();
      });
    });
  });

  apply();
  if (typeof window.registerPanelListener === "function") {
    window.registerPanelListener(function (s) {
      state = s || state;
      apply();
    });
  }
})();
