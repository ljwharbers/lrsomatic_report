// Genome circos driver labels: the tier 1-2 key findings of the ticked gene lists, written in the
// ring outside the chromosome names of the static circos SVG (R/circos.R).
//
// Inlined at the end of templates/per_sample.qmd after key_findings.js, whose kfRowMatches() is
// the chip test here too, so the labels and the overview table always list the same findings.
// window.CIRCOS_DRIVERS is written by templates/sections/_circos.qmd:
//   geom     {width, height, cx, cy, r_anchor, r_label, r_max} in SVG user units of the plate
//   anchors  [[gene, theta, x, y], ...]   theta: degrees, counter-clockwise from 3 o'clock
//   rows     [[label, group, genes, panels, tier, type, alteration, rescued], ...]  one per
//            finding line and member gene (circos_driver_rows() in R/key_findings.R); a
//            whitelist-rescued line (rescued 1) is labelled only while window.SV_SHOW_RESCUED
// The overlay <svg class="circos-drivers"> shares the plate image's viewBox, so the anchors R
// read off the open device land on the image exactly. Rebuilt from scratch on every listener
// call (registerPanelListener() fires more than once per change).
//
// The pure half (circosDriverLabels, placeLabels) is exported for node
// (tests/js/test-circos-drivers.js); the DOM half only runs when a document exists.
(function () {
  "use strict";

  // Minimum angular distance between two labels, degrees (label font ~7.5 units at r ~ 205)
  const MIN_GAP_DEG = 2.4;

  function split(x) { return String(x || "").split(",").filter(function (s) { return s.length > 0; }); }

  // Labels to show for a panel state: one per gene with an anchor, best tier first, carrying the
  // alterations of every matching line. rows/anchors as in window.CIRCOS_DRIVERS. A gene whose
  // every shown line is rescued is flagged `rescued`.
  function circosDriverLabels(rows, anchors, state, matches, showRescued) {
    const at = {};
    (anchors || []).forEach(function (a) { at[a[0]] = a; });
    const byGene = {}, order = [];
    (rows || []).forEach(function (r) {
      const tier = Number(r[4]);
      if (!(tier <= 2)) return;
      const rescued = Number(r[7]) === 1;
      if (rescued && !showRescued) return;
      const row = { gene: r[1], genes: split(r[2]).length ? split(r[2]) : [r[1]],
                    panels: split(r[3]), tier: tier };
      if (!matches(row, state)) return;
      const a = at[r[0]];
      if (!a) return;
      let g = byGene[r[0]];
      if (!g) {
        g = byGene[r[0]] = { gene: r[0], tier: tier, theta: Number(a[1]),
                             x: Number(a[2]), y: Number(a[3]), items: [], rescued: true };
        order.push(r[0]);
      }
      if (tier < g.tier) g.tier = tier;
      if (!rescued) g.rescued = false;
      const item = r[5] + ": " + r[6] + (rescued ? " [whitelist-rescued]" : "");
      if (g.items.indexOf(item) < 0) g.items.push(item);
    });
    return order.map(function (k) { return byGene[k]; });
  }

  // Clockwise degrees from 12 o'clock, [0, 360): the plot starts at the top (start.degree 90)
  // and its wide gap sits just left of it, so the seam at 0/360 falls in empty space
  function clockwise(theta) { return ((90 - theta) % 360 + 360) % 360; }

  // Spread labels so no two are closer than `gap` degrees while each cluster stays centred on
  // its genes: overlapping runs merge into one group laid out evenly around the mean of its
  // anchors, repeated until nothing overlaps. Sets `.labelTheta` on every item; returns them
  // in clockwise order.
  function placeLabels(items, gap) {
    gap = gap || MIN_GAP_DEG;
    const xs = items.map(function (it) { return { it: it, a: clockwise(it.theta) }; });
    xs.sort(function (p, q) { return p.a - q.a || (p.it.gene < q.it.gene ? -1 : 1); });
    let groups = xs.map(function (x) { return { m: [x], start: x.a }; });
    function layout(g) {
      const mean = g.m.reduce(function (s, x) { return s + x.a; }, 0) / g.m.length;
      const span = (g.m.length - 1) * gap;
      g.start = Math.min(Math.max(mean - span / 2, 0), Math.max(360 - span, 0));
    }
    let changed = true;
    while (changed && groups.length > 1) {
      changed = false;
      const next = [groups[0]];
      for (let i = 1; i < groups.length; i++) {
        const prev = next[next.length - 1], g = groups[i];
        if (g.start < prev.start + (prev.m.length - 1) * gap + gap - 1e-9) {
          prev.m = prev.m.concat(g.m);
          layout(prev);
          changed = true;
        } else next.push(g);
      }
      groups = next;
    }
    const out = [];
    groups.forEach(function (g) {
      g.m.forEach(function (x, k) {
        x.it.labelTheta = 90 - (g.start + k * gap);
        out.push(x.it);
      });
    });
    return out;
  }

  const api = { circosDriverLabels: circosDriverLabels, placeLabels: placeLabels,
                MIN_GAP_DEG: MIN_GAP_DEG };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  if (typeof window !== "undefined") window.circosDrivers = api;
  if (typeof document === "undefined") return;

  const D = window.CIRCOS_DRIVERS;
  const host = document.querySelector("svg.circos-drivers");
  const hint = document.querySelector(".circos-drivers-hint");
  if (!D || !host || !window.keyFindings) return;
  const G = D.geom;
  const NS = "http://www.w3.org/2000/svg";
  host.setAttribute("viewBox", "0 0 " + G.width + " " + G.height);

  function el(name, attrs) {
    const n = document.createElementNS(NS, name);
    Object.keys(attrs).forEach(function (k) { n.setAttribute(k, attrs[k]); });
    return n;
  }
  function pt(deg, r) {
    const t = deg * Math.PI / 180;
    return [G.cx + r * Math.cos(t), G.cy - r * Math.sin(t)];
  }
  function xy(p) { return p[0].toFixed(2) + "," + p[1].toFixed(2); }

  function draw(state) {
    while (host.firstChild) host.removeChild(host.firstChild);
    const labels = placeLabels(circosDriverLabels(D.rows, D.anchors, state,
                                                  window.keyFindings.kfRowMatches,
                                                  !!window.SV_SHOW_RESCUED));
    labels.forEach(function (g) {
      const cls = "tier-" + g.tier + (g.rescued ? " is-rescued" : "");
      const p0 = [g.x, g.y], p1 = pt(g.theta, G.r_anchor + 5), p2 = pt(g.labelTheta, G.r_label - 2);
      host.appendChild(el("polyline", { points: xy(p0) + " " + xy(p1) + " " + xy(p2),
                                        class: "circos-driver-line " + cls }));
      // Read outward on the right half, inward on the left, like circlize's niceFacing
      const flip = Math.cos(g.labelTheta * Math.PI / 180) < 0;
      const lp = pt(g.labelTheta, G.r_label);
      const rot = flip ? (180 - g.labelTheta) : -g.labelTheta;
      const t = el("text", {
        x: lp[0].toFixed(2), y: lp[1].toFixed(2),
        class: "circos-driver-label " + cls,
        "text-anchor": flip ? "end" : "start", "dominant-baseline": "middle",
        transform: "rotate(" + rot.toFixed(2) + " " + xy(lp).replace(",", " ") + ")",
        "data-gene": g.gene
      });
      t.textContent = g.gene;
      const title = el("title", {});
      title.textContent = g.gene + (g.rescued ? " (whitelist-rescued, tier " : " (tier ") +
                          g.tier + "): " + g.items.join("; ");
      t.appendChild(title);
      host.appendChild(t);
    });
    host.setAttribute("data-n-labels", String(labels.length));
    if (hint) {
      hint.textContent = state && !state.all
        ? (labels.length ? "" : "No tier 1–2 finding in the ticked gene list(s).")
        : "Tick a gene list to label its tier 1–2 findings.";
      hint.hidden = hint.textContent.length === 0;
    }
  }

  let last = { panels: [], customGenes: [], all: true };
  draw(last);
  if (typeof window.registerPanelListener === "function")
    window.registerPanelListener(function (s) { if (s) { last = s; draw(s); } });
  document.addEventListener("sv-rescue-change", function () { draw(last); });
})();
