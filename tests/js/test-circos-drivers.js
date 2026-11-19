// Genome circos driver labels (assets/js/circos_drivers.js): which findings get a label for a
// chip state, and the angular de-overlap.
//
// Plain `node`, no dependencies, like test-key-findings.js.
//
//   node tests/js/test-circos-drivers.js

const path = require("path");
const kf = require(path.join(__dirname, "..", "..", "assets", "js", "key_findings.js"));
const cd = require(path.join(__dirname, "..", "..", "assets", "js", "circos_drivers.js"));

let failed = 0;
function check(name, got, want) {
  const g = JSON.stringify(got), w = JSON.stringify(want);
  if (g !== w) {
    failed++;
    console.error("FAIL " + name + ": got " + g + ", want " + w);
  }
}

// [label, group, genes, panels, tier, type, alteration]
const rows = [
  ["EWSR1", "EWSR1–ATF1", "EWSR1,ATF1", "sarcoma", 1, "SV", "translocation"],
  ["ATF1",  "EWSR1–ATF1", "EWSR1,ATF1", "sarcoma", 1, "SV", "translocation"],
  ["CDKN2A", "CDKN2A", "CDKN2A", "sarcoma,lymphoid", 2, "SV", "deletion"],
  ["CDKN2A", "CDKN2A", "CDKN2A", "sarcoma,lymphoid", 1, "CNA", "homozygous deletion"],
  ["KMT2D", "KMT2D", "KMT2D", "lymphoid", 2, "small variant", "p.R5432*"],
  ["MYC",   "MYC",   "MYC",   "lymphoid", 3, "SV", "near MYC"],
  ["TP53",  "TP53",  "TP53",  "",         1, "small variant", "p.R175H"],
  ["NOPOS", "NOPOS", "NOPOS", "sarcoma",  1, "CNA", "amplification"]
];
// [gene, theta, x, y]
const anchors = [["EWSR1", -100, 1, 1], ["ATF1", 80, 2, 2], ["CDKN2A", 10, 3, 3],
                 ["KMT2D", 40, 4, 4], ["MYC", 60, 5, 5], ["TP53", 0, 6, 6]];
const st = (panels, custom, all) => ({ panels: panels, customGenes: custom || [], all: !!all });
const genes = (s) => cd.circosDriverLabels(rows, anchors, s, kf.kfRowMatches).map((g) => g.gene);

check("nothing ticked labels nothing", genes(st([], [], true)), []);
check("sarcoma: both partners of the joined row, best tier, no unanchored gene",
      genes(st(["sarcoma"])), ["EWSR1", "ATF1", "CDKN2A"]);
check("tier 3 is never labelled", genes(st(["lymphoid"])), ["CDKN2A", "KMT2D"]);
check("custom list reaches non-panel rows", genes(st([], ["TP53"])), ["TP53"]);
check("custom list matches the group's members", genes(st([], ["ATF1"])), ["EWSR1", "ATF1"]);

// Whitelist-rescued rows ([7] = 1) label only on request; a gene with a PASS line is not rescued
const rrows = rows.concat([
  ["NF1",   "NF1",   "NF1",   "sarcoma", 1, "SV", "translocation", 1],
  ["KMT2D", "KMT2D", "KMT2D", "lymphoid", 1, "SV", "deletion", 1]
]);
const ranchors = anchors.concat([["NF1", 20, 7, 7]]);
const rlab = (show) => cd.circosDriverLabels(rrows, ranchors, st(["sarcoma", "lymphoid"]),
                                             kf.kfRowMatches, show);
check("rescued: hidden by default", rlab(false).map((g) => g.gene).indexOf("NF1"), -1);
check("rescued: shown on request", rlab(true).map((g) => g.gene).indexOf("NF1") >= 0, true);
check("rescued: rescued-only gene flagged", rlab(true).filter((g) => g.gene === "NF1")[0].rescued, true);
const kmt = rlab(true).filter((g) => g.gene === "KMT2D")[0];
check("rescued: a PASS line keeps the gene unflagged, best tier wins", [kmt.rescued, kmt.tier], [false, 1]);
check("rescued: off, KMT2D keeps its PASS tier", rlab(false).filter((g) => g.gene === "KMT2D")[0].tier, 2);

const cdk = cd.circosDriverLabels(rows, anchors, st(["sarcoma"]), kf.kfRowMatches)
  .filter((g) => g.gene === "CDKN2A")[0];
check("one label per gene takes the best tier", cdk.tier, 1);
check("and lists every matching alteration", cdk.items,
      ["SV: deletion", "CNA: homozygous deletion"]);

// ---- placement -----------------------------------------------------------------------------
const cw = (theta) => ((90 - theta) % 360 + 360) % 360;
function gaps(placed) {
  const a = placed.map((p) => cw(p.labelTheta));
  const out = [];
  for (let i = 1; i < a.length; i++) out.push(a[i] - a[i - 1]);
  return out;
}
const gap = cd.MIN_GAP_DEG;

let placed = cd.placeLabels([{ gene: "A", theta: 30 }, { gene: "B", theta: -60 }]);
check("separate labels stay on their genes", placed.map((p) => p.labelTheta), [30, -60]);

placed = cd.placeLabels([{ gene: "A", theta: 30 }, { gene: "B", theta: 30.5 },
                         { gene: "C", theta: 29.5 }]);
check("a cluster is spread to the minimum gap",
      gaps(placed).every((d) => d >= gap - 1e-9), true);
const mean = placed.reduce((s, p) => s + cw(p.labelTheta), 0) / placed.length;
check("and stays centred on its genes", Math.abs(mean - cw(30)) < 1e-9, true);
check("clockwise order follows the genes", placed.map((p) => p.gene), ["B", "A", "C"]);

// Chains of clusters that only collide after a merge must merge too
const many = [];
for (let i = 0; i < 20; i++) many.push({ gene: "G" + i, theta: 45 - i * 0.3 });
many.push({ gene: "FAR", theta: 45 - 8 });
placed = cd.placeLabels(many);
check("chained clusters never overlap", gaps(placed).every((d) => d >= gap - 1e-9), true);

placed = cd.placeLabels([{ gene: "TOP", theta: 90.2 }, { gene: "TOP2", theta: 90 }]);
check("labels at the top seam stay within one turn",
      placed.every((p) => cw(p.labelTheta) >= 0 && cw(p.labelTheta) < 360), true);

if (failed) { console.error(failed + " check(s) failed"); process.exit(1); }
console.log("test-circos-drivers.js: all checks passed");
