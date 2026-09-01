// Row-matching tests for the tickbox column filters (assets/js/facet_filter.js).
//
// Dependency-free on purpose: plain `node`, no jsdom, no test framework — the repo has no
// JS toolchain and adding one would be a new dependency in recipe/meta.yaml's sibling
// environment.yml for no benefit. This covers the part of the widget that fails *silently*:
// the ext.search predicate. A wrong token split does not error, it just quietly shows the
// wrong rows. Everything DOM-shaped (menu placement, clipping, keyboard) still needs the
// browser checklist in the README/CLAUDE.md.
//
//   node tests/js/test-facet-predicate.js
//
// The stub below implements only what facet_filter.js touches at load time, plus
// ext.search.push, which is how the predicate is captured.

const fs = require("fs");
const path = require("path");
const vm = require("vm");

const SRC = path.join(__dirname, "..", "..", "assets", "js", "facet_filter.js");

const predicates = [];
const noop = () => {};
// A jQuery-shaped stand-in: every call returns a chainable object, so the module's
// document/window bindings and helper definitions run without a DOM.
function chain() {
  const obj = new Proxy(function () { return chain(); }, {
    get(_t, prop) {
      if (prop === "length") return 0;
      if (prop === "0") return undefined;
      return () => chain();
    },
    apply() { return chain(); }
  });
  return obj;
}
const $ = function () { return chain(); };
$.fn = { dataTable: { ext: { search: predicates }, isDataTable: () => false } };
$.trim = (s) => String(s).trim();

const sandbox = {
  $, jQuery: $, console,
  setTimeout, clearTimeout,
  // window.addEventListener is for the resize listener that re-places an open menu; the
  // document one is for the capturing scroll listener beside it.
  addEventListener: noop,
  document: { getElementById: () => null, addEventListener: noop },
};
sandbox.window = sandbox;

vm.createContext(sandbox);
vm.runInContext(fs.readFileSync(SRC, "utf8"), sandbox, { filename: SRC });

if (predicates.length !== 1) {
  console.error("expected exactly one ext.search predicate, got " + predicates.length);
  process.exit(1);
}
const predicate = predicates[0];
const STATE = sandbox.window.facetFilterState;

// Two tables, keyed by node identity exactly as the report keys them.
const snvNode = { table: "snv" };
const svNode  = { table: "sv" };
sandbox.window.snvTableElem = snvNode;
sandbox.window.svTableElem  = svNode;
sandbox.window.SNV_COLS   = { symbol: 0, consequence: 1, impact: 2, callers: 3 };
sandbox.window.SNV_FACETS = {
  consequence: { sep: ",", values: [] },
  impact:      { sep: null, values: [] },
  callers:     { sep: ",", values: [] }
};
sandbox.window.SV_COLS   = { id: 0, svclass: 1 };
sandbox.window.SV_FACETS = { svclass: { sep: null, values: [] } };

let failed = 0;
function check(name, got, want) {
  if (got !== want) {
    failed++;
    console.error("FAIL " + name + ": got " + got + ", want " + want);
  }
}
function tick(key, col, tokens) {
  STATE[key][col] = new Set(tokens);
}
function clear() {
  Object.keys(STATE).forEach((k) => Object.keys(STATE[k]).forEach((c) => delete STATE[k][c]));
}
// A row as DataTables hands it to a predicate: the raw, typed row array.
const snvRow = (consequence, impact, callers) =>
  ["TP53", consequence, impact, callers];
const run = (node, row) => predicate({ nTable: node }, row.map(String), 0, row);

// Nothing ticked: every row passes, and an unrelated DT table is never touched.
clear();
check("no selection passes", run(snvNode, snvRow("intron_variant", "MODIFIER", "")), true);
check("foreign table passes", predicate({ nTable: {} }, [], 0, ["x"]), true);

// A single ticked value on a whole-cell column.
clear();
tick("snv", "impact", ["HIGH"]);
check("impact HIGH matches",  run(snvNode, snvRow("stop_gained", "HIGH", "")), true);
check("impact LOW rejected",  run(snvNode, snvRow("synonymous_variant", "LOW", "")), false);

// OR within a column.
clear();
tick("snv", "impact", ["HIGH", "LOW"]);
check("impact OR first",  run(snvNode, snvRow("stop_gained", "HIGH", "")), true);
check("impact OR second", run(snvNode, snvRow("x", "LOW", "")), true);
check("impact OR other",  run(snvNode, snvRow("x", "MODERATE", "")), false);

// The multi-token contract: a ticked term matches a cell holding several.
clear();
tick("snv", "consequence", ["missense_variant"]);
check("multi-token cell matches on one term",
      run(snvNode, snvRow("missense_variant,splice_region_variant", "MODERATE", "")), true);
check("multi-token cell, other term ticked",
      run(snvNode, snvRow("intron_variant,non_coding_transcript_variant", "MODIFIER", "")),
      false);
check("whitespace around a token is ignored",
      run(snvNode, snvRow("splice_region_variant, missense_variant", "MODERATE", "")), true);

// A shared prefix must not match: this is substring matching's failure mode.
clear();
tick("snv", "consequence", ["splice_acceptor_variant"]);
check("prefix does not match",
      run(snvNode, snvRow("splice_acceptor_variant_x", "HIGH", "")), false);
check("exact token matches",
      run(snvNode, snvRow("splice_acceptor_variant", "HIGH", "")), true);

// A whole-cell column must NOT split on commas.
clear();
tick("snv", "impact", ["HIGH"]);
check("whole-cell column is not split", run(snvNode, snvRow("x", "HIGH,LOW", "")), false);

// The "(none)" bucket: null, undefined, "" and whitespace all land in it.
clear();
tick("snv", "impact", [null]);
check("none matches empty string",  run(snvNode, snvRow("x", "", "")), true);
check("none matches whitespace",    run(snvNode, snvRow("x", "   ", "")), true);
check("none matches null",          run(snvNode, snvRow("x", null, "")), true);
check("none matches undefined",     run(snvNode, snvRow("x", undefined, "")), true);
check("none rejects a real value",  run(snvNode, snvRow("x", "HIGH", "")), false);
clear();
tick("snv", "consequence", [null]);
check("none on a split column matches a comma-only cell",
      run(snvNode, snvRow(",,", "MODIFIER", "")), true);

// AND across columns.
clear();
tick("snv", "impact", ["HIGH"]);
tick("snv", "consequence", ["stop_gained"]);
check("both columns satisfied", run(snvNode, snvRow("stop_gained", "HIGH", "")), true);
check("one column fails",       run(snvNode, snvRow("stop_gained", "LOW", "")), false);
check("other column fails",     run(snvNode, snvRow("frameshift_variant", "HIGH", "")), false);

// Per-table isolation: the SV table's ticks must not filter the SNV table.
clear();
tick("sv", "svclass", ["translocation"]);
check("SV ticks leave SNV alone", run(snvNode, snvRow("intron_variant", "MODIFIER", "")), true);
check("SV row matches",  run(svNode, ["SV1", "translocation"]), true);
check("SV row rejected", run(svNode, ["SV2", "DEL"]), false);

// A column that vanished from *_COLS must not empty the table — the loud failure for that
// is the attach-time console.warn, not a silently blank table.
clear();
tick("snv", "nosuchcolumn", ["x"]);
check("unresolvable column is skipped",
      run(snvNode, snvRow("intron_variant", "MODIFIER", "")), true);

// ---- live counts ----------------------------------------------------------
// The number beside each value is accumulated inside the same predicate pass that does the
// matching, over the rows passing every *other* filter but ignoring the column's own ticks.
// It fails as silently as the matching does — a wrong number just reads as a wrong number —
// and it is the half that cannot be checked by eye on a 167k-row table.
//
// DataTables passes the per-predicate loop counter as the 5th argument, restarting at 0 for
// each predicate; the counts reset on that, so a scenario drives its first row with j = 0.
const runAt = (node, row, j) => predicate({ nTable: node }, row.map(String), 0, row, j);
const countOf = (key, col, token) => {
  const pairs = sandbox.window.facetFilterCounts()[key][col] || [];
  const hit = pairs.filter((p) => p[0] === token);
  return hit.length ? hit[0][1] : 0;
};
// Three rows exercising every shape at once: a whole-cell column (impact), a split column
// with a multi-token cell (consequence) and a split column with an empty cell (callers).
const COUNT_ROWS = [
  snvRow("missense_variant", "MODERATE", "clairs"),
  snvRow("intron_variant", "MODIFIER", "clairs,deepsomatic"),
  snvRow("missense_variant,splice_region_variant", "MODERATE", "")
];
function countPass() {
  COUNT_ROWS.forEach((row, j) => runAt(snvNode, row, j));
}

// Unfiltered: every row counts everywhere, a split cell adds to each of its tokens, and an
// empty cell lands in the (none) bucket. `consequence` sums to 4 over 3 rows, which is the
// over-sum js_facet_defs() documents on the R side.
clear();
countPass();
check("count missense (no filter)",     countOf("snv", "consequence", "missense_variant"), 2);
check("count intron (no filter)",       countOf("snv", "consequence", "intron_variant"), 1);
check("count splice_region (no filter)", countOf("snv", "consequence", "splice_region_variant"), 1);
check("count MODERATE (no filter)",     countOf("snv", "impact", "MODERATE"), 2);
check("count MODIFIER (no filter)",     countOf("snv", "impact", "MODIFIER"), 1);
check("count clairs (no filter)",       countOf("snv", "callers", "clairs"), 2);
check("count deepsomatic (no filter)",  countOf("snv", "callers", "deepsomatic"), 1);
check("count (none) bucket",            countOf("snv", "callers", null), 1);

// One column ticked: every other column follows it, and the ticked column's own counts do
// not — that is what keeps the remaining values of a ticked column readable. The row that
// fails only `impact` still counts towards `impact`, and towards nothing else.
clear();
tick("snv", "impact", ["MODERATE"]);
countPass();
check("ticked column ignores its own ticks (kept)",   countOf("snv", "impact", "MODERATE"), 2);
check("ticked column ignores its own ticks (excluded)", countOf("snv", "impact", "MODIFIER"), 1);
check("other column follows the tick (match)",  countOf("snv", "consequence", "missense_variant"), 2);
check("other column follows the tick (miss)",   countOf("snv", "consequence", "intron_variant"), 0);
check("split column follows the tick",          countOf("snv", "callers", "clairs"), 1);
check("split column follows the tick (zeroed)", countOf("snv", "callers", "deepsomatic"), 0);

// Two columns ticked: a row failing both counts nowhere at all, so MODIFIER — which the
// single-tick case above still credited — now reads 0.
clear();
tick("snv", "impact", ["MODERATE"]);
tick("snv", "consequence", ["missense_variant"]);
countPass();
check("row failing two columns counts nowhere", countOf("snv", "impact", "MODIFIER"), 0);
check("rows failing neither still count",       countOf("snv", "impact", "MODERATE"), 2);
check("companion token of a passing row",       countOf("snv", "consequence", "splice_region_variant"), 1);
check("value excluded by the other column",     countOf("snv", "consequence", "intron_variant"), 0);
check("unticked column sees both ticks",        countOf("snv", "callers", "clairs"), 1);

// A token repeated inside one cell counts its row once, mirroring
// js_facet_defs()'s unique(d, by = c("row", "tok")).
clear();
runAt(snvNode, snvRow("missense_variant", "HIGH", "clairs,clairs"), 0);
check("repeated token counts its row once", countOf("snv", "callers", "clairs"), 1);

// Counts are per table, exactly as the ticks are.
clear();
runAt(svNode, ["SV1", "translocation"], 0);
runAt(svNode, ["SV2", "DEL"], 1);
check("SV counts are the SV table's", countOf("sv", "svclass", "translocation"), 1);
check("SNV counts untouched by an SV pass", countOf("snv", "impact", "HIGH"), 1);

if (failed > 0) {
  console.error(failed + " check(s) failed");
  process.exit(1);
}
console.log("facet predicate: all checks passed");
