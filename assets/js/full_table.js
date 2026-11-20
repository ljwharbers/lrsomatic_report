// Full small-variant table embedded in the report (--embed-full-table): the sidecar's gzip bytes, base64, in an inert <script type="application/octet-stream"> that nothing reads until a click. "Download" hands those bytes back as <sample>_small_variants.tsv.gz (no decompression, works from file:// offline); "Show all" gunzips them with DecompressionStream, rebuilds DT rows from window.SNV_COL_SPEC and swaps them into the small-variant table. Pure helpers are exported for tests/js/test-full-table.js
(function () {
  "use strict";

  function b64ToBytes(b64) {
    var bin = atob(String(b64).replace(/[^A-Za-z0-9+/=]/g, ""));
    var out = new Uint8Array(bin.length);
    for (var i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
    return out;
  }

  function canDecompress(g) {
    g = g || (typeof window !== "undefined" ? window : globalThis);
    return typeof g.DecompressionStream === "function" && typeof g.Response === "function" &&
           typeof g.Blob === "function";
  }

  function gunzipToText(bytes) {
    var stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
    return new Response(stream).text();
  }

  // data.table::fwrite() TSV: quote = "auto" quotes a field only when it holds a tab, newline or quote, doubling inner quotes
  function parseTsv(text) {
    var rows = [];
    if (text.indexOf('"') < 0) {
      var lines = text.split("\n");
      for (var i = 0; i < lines.length; i++) {
        var ln = lines[i];
        if (ln.charCodeAt(ln.length - 1) === 13) ln = ln.slice(0, -1);
        if (ln === "" && i === lines.length - 1) continue;
        rows.push(ln.split("\t"));
      }
    } else {
      var row = [], field = "", q = false, n = text.length;
      for (var k = 0; k < n; k++) {
        var c = text[k];
        if (q) {
          if (c === '"') {
            if (text[k + 1] === '"') { field += '"'; k++; } else q = false;
          } else field += c;
        } else if (c === '"' && field === "") {
          q = true;
        } else if (c === "\t") {
          row.push(field); field = "";
        } else if (c === "\n") {
          if (field.charCodeAt(field.length - 1) === 13) field = field.slice(0, -1);
          row.push(field); rows.push(row); row = []; field = "";
        } else field += c;
      }
      if (field !== "" || row.length) { row.push(field); rows.push(row); }
    }
    return { header: rows.length ? rows[0] : [], rows: rows.slice(1) };
  }

  // DT escapes character cells (htmltools::htmlEscape) before embedding them; rows added here must match
  function htmlEscape(s) {
    return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

  // R's round() (nmath/fround.c, R >= 4): of the two candidates, the nearer in double arithmetic, an equal distance going to even; Math.round/toFixed differ on ~1 in 60 VAFs
  function roundR(x, d) {
    if (!isFinite(x) || x === 0) return x;
    if (x < 0) return -roundR(-x, d);
    var p = Math.pow(10, d), x10 = p * x, i10 = Math.floor(x10);
    var xd = i10 / p, xu = Math.ceil(x10) / p;
    var du = xu - x, dd = x - xd;
    return (du < dd || (du === dd && i10 % 2 === 1)) ? xu : xd;
  }

  function convert(v, type, scale, digits) {
    if (v === undefined || v === "" || v === "NA") return null;
    if (type === "number") {
      var x = Number(v);
      if (!isFinite(x)) return null;
      if (scale !== 1) {
        x = x * scale;
        if (digits !== null && digits !== undefined) x = roundR(x, digits);
      }
      return x;
    }
    if (type === "logical") return v === "TRUE" ? true : v === "FALSE" ? false : null;
    return htmlEscape(v);
  }

  // Sidecar rows -> DT row arrays placed by `cols` (name -> index, window.SNV_COLS); spec = [[display, source, type, scale, digits], ...] from snv_column_spec()
  function toTableRows(parsed, spec, cols) {
    var width = 0;
    Object.keys(cols).forEach(function (nm) { width = Math.max(width, cols[nm] + 1); });
    var plan = [];
    spec.forEach(function (s) {
      var dst = cols[s[0]];
      if (dst === undefined) return;
      plan.push([dst, parsed.header.indexOf(s[1]), s[2], Number(s[3]), s[4]]);
    });
    var out = new Array(parsed.rows.length);
    for (var r = 0; r < parsed.rows.length; r++) {
      var src = parsed.rows[r], row = new Array(width);
      for (var j = 0; j < width; j++) row[j] = null;
      for (var p = 0; p < plan.length; p++) {
        var pl = plan[p];
        row[pl[0]] = pl[1] < 0 ? null : convert(src[pl[1]], pl[2], pl[3], pl[4]);
      }
      out[r] = row;
    }
    return out;
  }

  function numericRange(rows, idx) {
    var lo = Infinity, hi = -Infinity;
    for (var r = 0; r < rows.length; r++) {
      var v = rows[r][idx];
      if (typeof v === "number" && isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
    }
    return lo <= hi ? [lo, hi] : null;
  }

  function downloadBytes(bytes, filename, doc, urlApi) {
    var blob = new Blob([bytes], { type: "application/gzip" });
    var url = urlApi.createObjectURL(blob);
    var a = doc.createElement("a");
    a.href = url;
    a.download = filename;
    a.style.display = "none";
    doc.body.appendChild(a);
    a.click();
    doc.body.removeChild(a);
    setTimeout(function () { urlApi.revokeObjectURL(url); }, 60000);
    return blob;
  }

  // The `swap` for facetFilterReload(): facets switch first, because widening a range slider fires DT's `set` handler, which draws
  function makeSwap(dt, rows, win, widen) {
    return function () {
      if (win.SNV_FACETS_FULL) win.SNV_FACETS = win.SNV_FACETS_FULL;
      dt.clear();
      dt.rows.add(rows);
      if (widen) widen();
    };
  }

  var api = { b64ToBytes: b64ToBytes, roundR: roundR, canDecompress: canDecompress, gunzipToText: gunzipToText,
              parseTsv: parseTsv, toTableRows: toTableRows, numericRange: numericRange,
              downloadBytes: downloadBytes, makeSwap: makeSwap };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  if (typeof document === "undefined") return;

  // DT's range sliders were sized to the embedded rows; widen them the way DT's own updateFilters() does (setFilterProps, minus Shiny)
  function widenSliders(dt, rows, spec, cols) {
    var $tds = $(dt.table().header()).children("tr").last().children("td");
    var vis = dt.columns(":visible").indexes().toArray();
    spec.forEach(function (s) {
      if (s[2] !== "number") return;
      var pos = vis.indexOf(cols[s[0]]);
      if (pos < 0) return;
      var $td = $tds.eq(pos);
      var slider = $td.find(".noUi-target").eq(0);
      if (!slider.length || typeof slider.noUiSlider !== "function") return;
      var rg = numericRange(rows, cols[s[0]]);
      if (!rg) return;
      try {
        var dec = Math.max(0, +slider.data("scale") || 0), scale = Math.pow(10, dec);
        var lims = [Math.floor(rg[0] * scale), Math.ceil(rg[1] * scale)];
        if (lims[0] === lims[1]) return;
        var vals = slider.val().map(Number), old = slider.noUiSlider("options").range;
        var start = lims.slice();
        if (vals[0] != old.min) start[0] = Math.max(vals[0], lims[0]);
        if (vals[1] != old.max) start[1] = Math.min(vals[1], lims[1]);
        slider.noUiSlider({ start: start, range: { min: lims[0], max: lims[1] } }, true);
        var fmt = function (x) {
          var v = (x / scale).toFixed(dec);
          return $td.attr("data-type") === "integer" ? String(parseInt(v, 10)) : v;
        };
        var $spans = slider.parent().children("span");
        $spans.first().text(fmt(start[0]));
        $spans.last().text(fmt(start[1]));
      } catch (err) {
        console.warn("full_table: could not widen the " + s[0] + " range filter.", err);
      }
    });
  }

  function init() {
    var holder = document.getElementById("snv-full-table-data");
    if (!holder) return;
    var box = holder.parentNode;
    var dl = box.querySelector(".full-table__download");
    var show = box.querySelector(".full-table__show");
    var status = box.querySelector(".full-table__status");
    var filename = holder.getAttribute("data-filename") || "small_variants.tsv.gz";
    var nRows = Number(holder.getAttribute("data-rows")) || 0;
    var bytes = null;
    var say = function (msg) { if (status) status.textContent = msg; };
    var payload = function () {
      if (!bytes) bytes = b64ToBytes(holder.textContent);
      return bytes;
    };

    if (dl) dl.addEventListener("click", function () {
      try {
        downloadBytes(payload(), filename, document, window.URL);
      } catch (err) {
        say("The download failed: " + err.message);
      }
    });

    if (!show) return;
    if (!canDecompress()) {
      say("This browser cannot unpack the table in place; the download works.");
      return;
    }
    show.hidden = false;
    show.addEventListener("click", function () {
      var node = window.snvTableElem;
      var spec = window.SNV_COL_SPEC, cols = window.SNV_COLS;
      if (!node || !$.fn.dataTable.isDataTable(node) || !spec || !cols) {
        say("The table is not ready yet.");
        return;
      }
      show.disabled = true;
      say("Loading " + nRows.toLocaleString() + " variants…");
      gunzipToText(payload()).then(function (text) {
        var rows = toTableRows(parseTsv(text), spec, cols);
        var dt = $(node).DataTable();
        var swap = makeSwap(dt, rows, window, function () { widenSliders(dt, rows, spec, cols); });
        if (window.facetFilterReload) window.facetFilterReload("snv", swap);
        else { swap(); dt.draw(false); }
        var scope = box.parentNode && box.parentNode.querySelector(".full-table__scope");
        if (scope) scope.textContent = "Showing all " + rows.length.toLocaleString() +
                                       " small variants.";
        show.hidden = true;
        say("");
      }).catch(function (err) {
        show.disabled = false;
        say("Could not unpack the table (" + err.message + "); the download still works.");
      });
    });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
