// Row-visibility tests for the key-findings table (assets/js/key_findings.js).
//
// Plain `node`, no dependencies, like test-facet-predicate.js. The module exports its pure
// predicates when `module` exists and skips its DOM half when `document` does not, so it is
// required directly. Covered: panel sets, custom genes, the empty ("all") state, tier-3
// collapsing, and idempotence (the panel listener fires more than once per change).
//
//   node tests/js/test-key-findings.js

const path = require("path");
const kf = require(path.join(__dirname, "..", "..", "assets", "js", "key_findings.js"));

let failed = 0;
function check(name, got, want) {
  const g = JSON.stringify(got), w = JSON.stringify(want);
  if (g !== w) {
    failed++;
    console.error("FAIL " + name + ": got " + g + ", want " + w);
  }
}

const rows = [
  { gene: "EWSR1", panels: ["sarcoma", "lymphoid"], tier: 1 },
  { gene: "ATF1",  panels: ["sarcoma"],             tier: 1 },
  { gene: "KMT2D", panels: ["lymphoid"],            tier: 2 },
  { gene: "MYC",   panels: ["lymphoid"],            tier: 3 },
  { gene: "TP53",  panels: [""],                    tier: 1 }   // non-panel row
];
const st = (panels, custom, all) => ({ panels: panels, customGenes: custom || [], all: !!all });

// Empty state: nothing ticked shows nothing, whatever the rows
let s = kf.kfSummarise(rows, st([], [], true), false);
check("all: nothing visible", s.visible, [false, false, false, false, false]);
check("all: flagged", s.all, true);
check("all: not 'empty'", s.empty, false);
check("null state is all", kf.kfSummarise(rows, null, true).visible.some(Boolean), false);

// One panel
s = kf.kfSummarise(rows, st(["sarcoma"]), false);
check("sarcoma rows", s.visible, [true, true, false, false, false]);
check("sarcoma counts", [s.matched, s.shown, s.nLow], [2, 2, 0]);

// Union of two panels; tier 3 collapsed until showLow
s = kf.kfSummarise(rows, st(["sarcoma", "lymphoid"]), false);
check("union rows, tier 3 hidden", s.visible, [true, true, true, false, false]);
check("union nLow", s.nLow, 1);
s = kf.kfSummarise(rows, st(["sarcoma", "lymphoid"]), true);
check("union rows, tier 3 shown", s.visible, [true, true, true, true, false]);

// A panel with no findings is the "empty" state, not the "all" state
s = kf.kfSummarise(rows, st(["other"]), false);
check("unmatched panel: empty", [s.empty, s.all, s.matched], [true, false, 0]);

// Custom genes match by symbol, case-insensitively, including non-panel rows
s = kf.kfSummarise(rows, st([], ["tp53", "MYC"]), false);
check("custom: TP53 shown, MYC tier 3 hidden", s.visible, [false, false, false, false, true]);
check("custom: MYC counted as low", s.nLow, 1);
s = kf.kfSummarise(rows, st(["sarcoma"], ["KMT2D"]), false);
check("custom + panel union", s.visible, [true, true, true, false, false]);

// Panels arrive as arrays from the listener but a comma string must work too
check("string panels", kf.kfRowMatches({ gene: "X", panels: "a,sarcoma", tier: 1 }, st(["sarcoma"])), true);
// An empty panel key on the row must never match an empty key in the state
check("empty key never matches", kf.kfRowMatches({ gene: "X", panels: [""], tier: 1 }, st([""])), false);

// Idempotence: the same state twice gives the same answer, and nothing is mutated
const state = st(["lymphoid"], ["ATF1"]);
const a = kf.kfSummarise(rows, state, false);
const b = kf.kfSummarise(rows, state, false);
check("idempotent", a, b);
check("state untouched", state, st(["lymphoid"], ["ATF1"]));

// Print label
check("label all", kf.kfListLabel(st([], [], true), {}), "none (no gene list ticked)");
check("label panels", kf.kfListLabel(st(["lymphoid", "sarcoma"]), { lymphoid: "Lymphoid" }), "Lymphoid, sarcoma");
check("label custom", kf.kfListLabel(st([], ["A", "B"]), {}), "Custom (2 genes)");

if (failed) {
  console.error(failed + " check(s) failed");
  process.exit(1);
}
console.log("key findings: all checks passed");
