// Clinical overview: the key-findings table follows the gene-list chips.
//
// Inlined at the end of templates/per_sample.qmd by the key-findings-js chunk (same route as
// facet_filter.js), i.e. after every section's markup and after the panel-js-data globals.
// The markup is templates/sections/_overview.qmd: one <tbody class="kf-group"> per gene (or
// per rearrangement joining two listed genes, "EWSR1–ATF1"), holding one <tr class="kf-row">
// per alteration of that gene:
//   tbody  data-gene   group key            data-genes   comma-joined member symbols
//   tr     data-panels comma-joined panel keys ("" = no panel)
//          data-tier   1 | 2 | 3            data-type    small variant | SV | CNA
// Rows are computed in R for EVERY embedded panel (R/key_findings.R); this script only shows
// and hides them. A line is shown when its panels intersect the ticked ones, or one of its
// group's genes is in a typed Custom… list; a group is shown when any of its lines is. The
// first shown line of a group is its head (R sorts lines by tier, so the head carries the
// group's best shown tier) and is the only one that prints tier, gene and lists; groups are
// re-ordered by best shown tier, R's order breaking ties. Nothing ticked ("all") shows no
// findings at all, but a notice with per-list totals that tick the chip when clicked: a
// findings table without a gene list would be every HIGH-impact variant in the genome.
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

  // Does one finding line match the panel state? row = {gene, genes: [...], panels: [...], tier}
  function kfRowMatches(row, state) {
    if (!state || state.all) return false;
    const panels = asList(state.panels);
    const rowPanels = asList(row.panels).filter(function (p) { return p.length > 0; });
    for (let i = 0; i < rowPanels.length; i++) {
      if (panels.indexOf(rowPanels[i]) >= 0) return true;
    }
    const custom = asList(state.customGenes).map(function (g) { return String(g).toUpperCase(); });
    if (custom.length === 0) return false;
    const genes = asList(row.genes).length ? asList(row.genes) : [row.gene || ""];
    return genes.some(function (g) { return custom.indexOf(String(g).toUpperCase()) >= 0; });
  }

  // Visibility of every line and group plus the counts the widget prints. Tier-3 lines that
  // match stay hidden unless showLow; they are counted in nLow so the toggle can say how many.
  // row.group is the group index (defaults to one group per row).
  //   groups[k] = {rows, matched, visible, head (row index | -1), tier, panels (shown union)}
  //   order     = group indices, shown groups by best shown tier, then R order
  function kfSummarise(rows, state, showLow) {
    const all = !state || !!state.all;
    const visible = [];
    const groups = [], byKey = {};
    let matched = 0, nLow = 0, shown = 0;
    rows.forEach(function (row, i) {
      const m = kfRowMatches(row, state);
      const low = Number(row.tier) >= 3;
      if (m) matched++;
      if (m && low) nLow++;
      const v = m && (!low || !!showLow);
      if (v) shown++;
      visible.push(v);

      const key = row.group === undefined ? "r" + i : String(row.group);
      if (!(key in byKey)) {
        byKey[key] = groups.length;
        groups.push({ rows: [], matched: false, visible: false, head: -1,
                      tier: Infinity, panels: [] });
      }
      const g = groups[byKey[key]];
      g.rows.push(i);
      if (m) g.matched = true;
      if (v) {
        g.visible = true;
        const t = Number(row.tier);
        if (g.head < 0 || t < g.tier) { g.head = i; g.tier = t; }
        asList(row.panels).forEach(function (p) {
          if (p.length > 0 && g.panels.indexOf(p) < 0) g.panels.push(p);
        });
      }
    });
    const order = groups.map(function (g, k) { return k; });
    order.sort(function (a, b) {
      const ta = groups[a].visible ? groups[a].tier : Infinity;
      const tb = groups[b].visible ? groups[b].tier : Infinity;
      if (ta !== tb) return ta < tb ? -1 : 1;
      return a - b;
    });
    const nGroups = groups.filter(function (g) { return g.matched; }).length;
    const nGroupsShown = groups.filter(function (g) { return g.visible; }).length;
    return { visible: visible, matched: matched, nLow: nLow, shown: shown, all: all,
             empty: !all && matched === 0, groups: groups, order: order,
             nGroups: nGroups, nGroupsShown: nGroupsShown };
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
  const table = root.querySelector("table.kf-table");
  const tbodies = Array.prototype.slice.call(root.querySelectorAll("tbody.kf-group"));
  const trs = [], rows = [];
  tbodies.forEach(function (tb, k) {
    const genes = (tb.getAttribute("data-genes") || tb.getAttribute("data-gene") || "").split(",");
    Array.prototype.forEach.call(tb.querySelectorAll("tr.kf-row"), function (tr) {
      trs.push(tr);
      rows.push({ gene: tr.getAttribute("data-gene") || "", genes: genes, group: k,
                  panels: (tr.getAttribute("data-panels") || "").split(","),
                  tier: Number(tr.getAttribute("data-tier")) });
    });
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
    trs.forEach(function (tr, i) {
      tr.hidden = !s.visible[i];
      const head = s.groups[rows[i].group].head === i;
      tr.classList.toggle("kf-row--head", head);
      tr.classList.toggle("kf-row--cont", s.visible[i] && !head);
    });
    s.groups.forEach(function (g, k) {
      tbodies[k].hidden = !g.visible;
      if (g.head < 0) return;
      // The head prints the lists of every shown line of its group
      const cell = trs[g.head].querySelector("td.kf-lists");
      if (cell) {
        cell.textContent = g.panels.length
          ? g.panels.map(function (p) { return labels[p] || p; }).join(", ") : "—";
      }
    });
    // Re-order only when it changed: appendChild moves the node, so this is idempotent
    if (table) {
      const current = Array.prototype.slice.call(table.querySelectorAll("tbody.kf-group"));
      const wanted = s.order.map(function (k) { return tbodies[k]; });
      if (current.some(function (tb, i) { return tb !== wanted[i]; })) {
        wanted.forEach(function (tb) { table.appendChild(tb); });
      }
    }
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
        (s.nGroupsShown + " of " + s.nGroups + " gene" + (s.nGroups === 1 ? "" : "s") + " shown");
    }
    if (printLine) printLine.textContent = "Gene list(s): " + kfListLabel(state, labels);
  }

  if (toggle) toggle.addEventListener("click", function () { showLow = !showLow; apply(); });

  // "Lymphoid · 4 genes" hints in the empty state tick the chip itself, through its own change
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
