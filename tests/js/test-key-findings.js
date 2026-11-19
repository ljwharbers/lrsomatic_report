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

// ---- gene groups ---------------------------------------------------------------------------
// R order: group 0 = KRAS (tier 2 only), group 1 = CDKN2A (CN tier 1 in sarcoma, SV tier 3
// in lymphoid), group 2 = EWSR1–ATF1 joined row (genes from the tbody), group 3 = MYC tier 1
const grows = [
  { gene: "KRAS",       genes: ["KRAS"],          group: 0, panels: ["lymphoid"],            tier: 2 },
  { gene: "CDKN2A",     genes: ["CDKN2A"],        group: 1, panels: ["sarcoma"],             tier: 1 },
  { gene: "CDKN2A",     genes: ["CDKN2A"],        group: 1, panels: ["lymphoid"],            tier: 3 },
  { gene: "EWSR1–ATF1", genes: ["EWSR1", "ATF1"], group: 2, panels: ["sarcoma", "lymphoid"], tier: 1 },
  { gene: "MYC",        genes: ["MYC"],           group: 3, panels: ["lymphoid"],            tier: 1 }
];
s = kf.kfSummarise(grows, st(["lymphoid"]), false);
check("groups: KRAS, joined, MYC shown; CDKN2A's tier-3 line hidden",
      s.groups.map(function (g) { return g.visible; }), [true, false, true, true]);
check("groups: shown groups by best tier, then R order", s.order, [2, 3, 0, 1]);
check("groups: counts", [s.nGroups, s.nGroupsShown], [4, 3]);
s = kf.kfSummarise(grows, st(["lymphoid", "sarcoma"]), true);
check("groups: head is the best-tier shown line", s.groups[1].head, 1);
check("groups: lists are the union of shown lines", s.groups[1].panels, ["sarcoma", "lymphoid"]);
s = kf.kfSummarise(grows, st(["lymphoid"]), true);
check("groups: a tier-3-only group sorts last", s.order, [2, 3, 0, 1]);
check("groups: head of CDKN2A under lymphoid is its tier-3 line", [s.groups[1].head, s.groups[1].tier], [2, 3]);
check("groups: hidden lines do not add lists", s.groups[1].panels, ["lymphoid"]);
// Custom… matches any member gene of a joined row
s = kf.kfSummarise(grows, st([], ["atf1"]), false);
check("custom: joined row matched by its second gene", s.visible, [false, false, false, true, false]);
check("custom: kfRowMatches on genes", kf.kfRowMatches(grows[3], st([], ["EWSR1"])), true);
// A row without genes falls back to gene
check("custom: no genes, gene used", kf.kfRowMatches({ gene: "TP53", panels: [""], tier: 1 }, st([], ["TP53"])), true);
// Rows without a group are one group each (old markup)
check("ungrouped rows: one group per row", kf.kfSummarise(rows, st(["sarcoma"]), false).groups.length, rows.length);
// Idempotent with groups too
const gs = st(["sarcoma"], ["MYC"]);
check("groups idempotent", kf.kfSummarise(grows, gs, true), kf.kfSummarise(grows, gs, true));

// ---- whitelist-rescued SVs -----------------------------------------------------------------
// Hidden until showRescued, whatever showLow says; counted apart from tier 3
const rrows = [
  { gene: "EWSR1", genes: ["EWSR1"], group: 0, panels: ["sarcoma"], tier: 1 },
  { gene: "EWSR1", genes: ["EWSR1"], group: 0, panels: ["sarcoma"], tier: 4, rescued: true },
  { gene: "NF1",   genes: ["NF1"],   group: 1, panels: ["sarcoma"], tier: 4, rescued: true },
  { gene: "KIT",   genes: ["KIT"],   group: 2, panels: ["sarcoma"], tier: 3 }
];
s = kf.kfSummarise(rrows, st(["sarcoma"]), true, false);
check("rescued: hidden by default, tier 3 follows showLow", s.visible, [true, false, false, true]);
check("rescued: counted apart", [s.nRescued, s.nLow], [2, 1]);
check("rescued: a rescued-only gene is not shown", s.groups[1].visible, false);
s = kf.kfSummarise(rrows, st(["sarcoma"]), false, true);
check("rescued: shown on request, tier 3 still collapsed", s.visible, [true, true, true, false]);
check("rescued: the gene head stays the PASS line", s.groups[0].head, 0);
check("rescued: rescued-only group sorts after tier 1", s.order, [0, 1, 2]);
check("rescued: unmatched lists count nothing", kf.kfSummarise(rrows, st(["lymphoid"]), true, true).nRescued, 0);

// Print label
check("label all", kf.kfListLabel(st([], [], true), {}), "none (no gene list ticked)");
check("label panels", kf.kfListLabel(st(["lymphoid", "sarcoma"]), { lymphoid: "Lymphoid" }), "Lymphoid, sarcoma");
check("label custom", kf.kfListLabel(st([], ["A", "B"]), {}), "Custom (2 genes)");

if (failed) {
  console.error(failed + " check(s) failed");
  process.exit(1);
}
console.log("key findings: all checks passed");
