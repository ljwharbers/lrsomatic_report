// Checkbox-dropdown ("tickbox") column filters for the small-variant and SV tables.
//
// DT gives every column a free-text search box (filter = "top"), which only helps a reader
// who already knows what the column contains — nobody remembers the ~30 VEP consequence
// terms. For the categorical columns this replaces that box with a dropdown listing the
// values actually present in this sample, with per-value row counts, so the filter doubles
// as a summary of the column.
//
// R publishes the values (window.SNV_FACETS / window.SV_FACETS, from js_facet_defs() in
// R/utils.R) and the column-name -> index maps (window.SNV_COLS / window.SV_COLS) that
// every by-name lookup here goes through. R's counts are the *unfiltered baseline* (and the
// value ordering); the number shown beside each value is recomputed on every filtering pass
// — see "live counts" below. Five things about the DOM and about DataTables this works in,
// each learned from internals rather than assumed:
//
//   1. The button goes in DT's own `div.form-group`, which is KEPT with its <input> only
//      hidden (DT holds a reference to that input, and the div is what DT's own
//      $td.children('div').first()/.last() lookups find). The MENU is not in the cell at
//      all: it is a child of <body>, position:fixed, placed each time from the button's
//      getBoundingClientRect(). Both halves of that are forced —
//        * the cell sits inside .dataTables_scrollHead and .dataTables_scrollBody, which
//          clip an absolutely-positioned menu. DT ships a workaround for its own filters
//          (flip those containers to overflow:visible on a "show" event, datatables.js
//          ~474-500) and it is unusable here: .dataTables_scrollBody is what *holds* the
//          horizontal scroll position, so making it visible discards it — on a wide table
//          the whole view snaps back to the first column and the menu you just opened is
//          off screen. We therefore never trigger DT's show/hide.
//        * a fixed-position menu left inside the cell would be duplicated, visibly, by the
//          sizing clone in note 2: that clone clips an absolute child in a zero-height
//          wrapper, but nothing clips a fixed one — and ticking a value causes a draw,
//          which is when the clone is rebuilt.
//   2. No `id` attributes anywhere in this markup, and no handlers delegated from
//      `document` onto widget classes. DataTables rebuilds a zero-height sizing clone of
//      the whole header in the scroll body on every draw, refilling each cell by innerHTML
//      after stripping ids — an id here would reappear duplicated once per draw. Hence
//      aria-expanded rather than aria-controls, handlers bound to the nodes we create, and
//      every query rooted at api.table().header(), never at `document`. The count nodes
//      cached for the live update are in that live header, which DT moves but never rebuilds.
//   3. The R-side column index is NOT the DOM position: DataTables detaches the filter
//      cells of visible:false columns. columns(':visible').indexes() translates it, and the
//      header label above the chosen cell is compared to the column name — a reorder that
//      broke the translation would otherwise attach a working dropdown to the wrong column.
//   4. Ticked state lives here, in a closure, not in the DOM and not in DT's per-column
//      search slots. That is what makes it survive a gene-panel redraw, Scroller's node
//      recycling and the per-draw header clone.
//   5. This predicate must be the LAST entry in $.fn.dataTable.ext.search. DataTables'
//      _fnFilterCustom runs the predicates in push order, each narrowing the previous
//      result, and only after the global search box and every per-column text box — so
//      being last is exactly what makes the rows we see "everything the reader has already
//      filtered down to", the gene-panel predicate (per_sample.qmd) included. attach()
//      re-pushes us to the end rather than relying on script order to stay as it is.
//
// ---- live counts -----------------------------------------------------------------------
// The count beside a value answers "what would I get if I ticked this", so it is measured
// over the rows passing every *other* filter — the gene panel, the text boxes and the other
// facet columns — while ignoring that column's own ticks. Ticking one consequence therefore
// leaves the other consequence counts readable (the standard facet-search convention);
// zeroed values grey out in place rather than moving, since the order stays R's.
//
// It is accumulated inside the matching pass, not in a scan of its own: a row that fails no
// column belongs to every column's count, a row that fails exactly one belongs to that
// column's only, and a row that fails two belongs to none. What makes that affordable on a
// 167k-row table is that the tokens are precomputed once, into a per-column CSR index of
// integer ids (starts/ids Int32Arrays keyed by DataTables row index) — splitting three
// columns of 167k rows on every pass is the ~150-500 ms cost that made the counts static in
// the first place. Splitting live is kept as the fallback for the pass that runs before
// attach() has built the index, and for tests/js/test-facet-predicate.js, which has no
// DataTables instance to build one from.
(function () {
  "use strict";

  // consequence is ~30 terms; impact is 4. Below this a find box is just clutter.
  var SEARCH_THRESHOLD = 12;
  // One draw over 167k rows runs two ext.search predicates per row, so coalesce rapid ticks.
  var DRAW_DEBOUNCE_MS = 120;
  var LABEL_MAX = 16;

  // Which tables get facets, keyed on the globals they publish from initComplete — the same
  // identity test the gene-panel predicate in per_sample.qmd uses, which excludes the QC,
  // ASCAT and phasing DT tables without any extra condition.
  var TABLES = [
    { key: "snv",
      elem:      function () { return window.snvTableElem; },
      cols:      function () { return window.SNV_COLS; },
      facets:    function () { return window.SNV_FACETS; },
      requested: function () { return window.SNV_FACET_COLS; } },
    { key: "sv",
      elem:      function () { return window.svTableElem; },
      cols:      function () { return window.SV_COLS; },
      facets:    function () { return window.SV_FACETS; },
      requested: function () { return window.SV_FACET_COLS; } }
  ];

  // key -> { colName -> Set of ticked tokens }. `null` in a Set is the "(none)" bucket.
  // A column absent from a table's object is unconstrained.
  var STATE    = { snv: {}, sv: {} };
  var ATTACHED = {};
  var WIDGETS  = [];      // every built widget, for the "clear all" control
  var PENDING  = {};      // per-table draw debounce timers
  var OPEN     = null;    // the one open menu

  // Counting machinery, all keyed by table key. FCOLS is the facet columns that resolved in
  // *_COLS, in payload order; everything else is parallel to it by column name.
  var FCOLS   = {};       // key -> [colName]
  var VOCAB   = {};       // key -> colName -> { strs: [token|null], map: Map, noneId: int }
  var CIDX    = {};       // key -> colName -> { starts: Int32Array, ids: Int32Array }
  var COUNTS  = {};       // key -> colName -> [rows per vocabulary id]
  var TOUCHED = {};       // key -> has the predicate run since the last flush?

  // Per-pass scratch, indexed by position in FCOLS. Reused rather than reallocated: the
  // predicate runs once per row, and only one table filters at a time.
  var SC_SRC = [], SC_LO = [], SC_HI = [], SC_TMP = [];

  function tableFor(node) {
    for (var i = 0; i < TABLES.length; i++) if (TABLES[i].elem() === node) return TABLES[i];
    return null;
  }

  function esc(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;")
                    .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }

  // ---- matching ------------------------------------------------------------
  // Mirrors how the cell was built, which is why the separator comes from R rather than
  // being guessed here: `consequence` is VEP's "&"-joined terms rewritten to commas
  // (R/parse_smallvariants.R) and `callers` is paste(sort(unique(caller)), collapse = ",").
  // NA and "" both collapse to the single null bucket the dropdown labels "(none)".
  function tokens(v, sep) {
    if (v === null || v === undefined) return [null];
    var s = String(v).trim();
    if (s === "") return [null];
    if (!sep) return [s];
    var out = s.split(sep).map(function (x) { return x.trim(); })
                          .filter(function (x) { return x.length > 0; });
    return out.length ? out : [null];
  }

  // ---- vocabulary ----------------------------------------------------------
  // Token strings are interned to integer ids so a count is an array increment rather than
  // a Map lookup. Seeded from the payload, in payload order, so a widget's option i and
  // vocabulary id i coincide; it still auto-extends, because the fallback split path can
  // meet a token js_facet_defs() did not list (it cannot happen with a matching payload,
  // but a silently mis-sized array would be worse than a slightly longer one).
  function addTok(v, tok) {
    var id = v.strs.length;
    v.strs.push(tok);
    if (tok === null) v.noneId = id; else v.map.set(tok, id);
    return id;
  }

  function idOf(v, tok) {
    if (tok === null || tok === undefined) {
      return v.noneId >= 0 ? v.noneId : addTok(v, null);
    }
    var id = v.map.get(tok);
    return id === undefined ? addTok(v, tok) : id;
  }

  function newVocab(def) {
    var v = { strs: [], map: new Map(), noneId: -1 };
    var vals = (def && def.values) || [];
    for (var i = 0; i < vals.length; i++) addTok(v, vals[i][0]);
    return v;
  }

  // Distinct ids of one cell, written into `out`; returns how many. Deduplicated per row so
  // "a,a" counts its row once, exactly as js_facet_defs()'s unique(d, by = c("row","tok")).
  function idsInto(voc, toks, out) {
    var n = 0;
    for (var i = 0; i < toks.length; i++) {
      var id = idOf(voc, toks[i]), dup = false;
      for (var k = 0; k < n; k++) if (out[k] === id) { dup = true; break; }
      if (!dup) out[n++] = id;
    }
    return n;
  }

  function resetCounts(key) {
    var cnt = COUNTS[key], voc = VOCAB[key];
    if (!cnt) return;
    Object.keys(cnt).forEach(function (nm) {
      var a = cnt[nm], n = voc[nm].strs.length;
      a.length = n;
      for (var i = 0; i < n; i++) a[i] = 0;
    });
  }

  // Set up the counting structures for a table the first time it is seen. Called from the
  // predicate as well as from attach(), so counts exist even before the widgets do.
  function ensureTable(t) {
    if (FCOLS[t.key]) return FCOLS[t.key].length > 0;
    var defs = t.facets(), cols = t.cols();
    if (!defs || !cols) return false;             // not published yet — retry on the next row
    var names = [], voc = {}, cnt = {};
    Object.keys(defs).forEach(function (nm) {
      if (cols[nm] === undefined) return;         // renamed column: attach() warns loudly
      names.push(nm);
      voc[nm] = newVocab(defs[nm]);
      cnt[nm] = [];
    });
    FCOLS[t.key]  = names;
    VOCAB[t.key]  = voc;
    COUNTS[t.key] = cnt;
    resetCounts(t.key);
    return names.length > 0;
  }

  function bump(arr, src, lo, hi) {
    for (var k = lo; k < hi; k++) {
      var id = src[k];
      arr[id] = (arr[id] || 0) + 1;
    }
  }

  // One tokenisation of every row of every facet column, into a CSR layout keyed by the
  // DataTables row index the predicate is handed. ~6 MB for the 167k-row SNV table, against
  // an array-of-arrays that would be an order of magnitude more. Public API only —
  // rows().indexes() and rows().data() are in the same order, and no row is ever added to
  // or removed from these tables, so the index stays valid for the session.
  function buildIndex(t, api) {
    var names = FCOLS[t.key];
    if (!names || !names.length) return;
    var defs = t.facets(), cols = t.cols();
    var idxs = api.rows().indexes().toArray();
    var data = api.rows().data();
    var maxRow = -1, k;
    for (k = 0; k < idxs.length; k++) if (idxs[k] > maxRow) maxRow = idxs[k];
    if (maxRow < 0) return;

    var byIdx = new Array(maxRow + 1);
    for (k = 0; k < idxs.length; k++) byIdx[idxs[k]] = data[k];

    var idx = {}, scratch = [];
    for (var c = 0; c < names.length; c++) {
      var nm = names[c], voc = VOCAB[t.key][nm], sep = defs[nm].sep, col = cols[nm];
      var starts = new Int32Array(maxRow + 2);
      var buf = new Int32Array(Math.max(16, (maxRow + 1) * 2));
      var total = 0;
      for (var r = 0; r <= maxRow; r++) {
        var row = byIdx[r];
        var n = (row === undefined) ? 0 : idsInto(voc, tokens(row[col], sep), scratch);
        if (total + n > buf.length) {
          var grown = new Int32Array(Math.max(buf.length * 2, total + n));
          grown.set(buf);
          buf = grown;
        }
        for (var q = 0; q < n; q++) buf[total + q] = scratch[q];
        total += n;
        starts[r + 1] = total;
      }
      idx[nm] = { starts: starts, ids: buf.subarray(0, total) };
    }
    CIDX[t.key] = idx;
  }

  // One predicate for both tables, keyed on the table node exactly as the gene-panel
  // predicate is. OR within a column, AND across columns. Both predicates are ext.search
  // entries, so DataTables ANDs them: a tick composes with the panel filter and with every
  // remaining per-column text box for free.
  //
  // Chosen over column().search() with a regex because the token rule is the load-bearing
  // part: as a regex it would be re-encoded in a second language with an escaping
  // obligation (svclass values carry spaces, callers carries "clairs-to"), and an
  // unescaped metacharacter yields a filter that matches MORE rows and says nothing.
  //
  // The exclude-own-column counts are accumulated here too, in the same visit — see the
  // header. `counter` is DataTables' per-predicate loop index, which restarts at 0 for each
  // predicate and so marks the first row of a pass.
  function facetPredicate(settings, searchData, index, rowData, counter) {
    var t = tableFor(settings.nTable);
    if (!t) return true;                                  // another DT table in the report
    if (!ensureTable(t)) return true;                     // no facet column on this table

    var key = t.key;
    if (counter === 0 || !TOUCHED[key]) resetCounts(key);
    TOUCHED[key] = true;

    var names = FCOLS[key], sel = STATE[key], voc = VOCAB[key], cnt = COUNTS[key];
    var cols = t.cols(), defs = t.facets(), idx = CIDX[key];
    // rowData is the original (typed) row, searchData the rendered strings. Always the
    // former: the SV table renders `locus`/`size` orthogonally, and mixing the two sources
    // per column is the trap variant_key() documents on the R side.
    var row = rowData || searchData;
    var fails = 0, failed = -1, i;

    for (i = 0; i < names.length; i++) {
      var nm = names[i], chosen = sel[nm];
      var hasSel = !!(chosen && chosen.size);
      // Once one column has failed, a column with nothing ticked can neither fail nor be
      // counted for this row — only the single column it failed can be.
      if (!hasSel && fails > 0) { SC_LO[i] = SC_HI[i] = 0; continue; }

      var src, lo, hi;
      var ci = idx && idx[nm];
      if (ci && index >= 0 && index + 1 < ci.starts.length) {
        src = ci.ids; lo = ci.starts[index]; hi = ci.starts[index + 1];
      } else {
        var out = SC_TMP[i] || (SC_TMP[i] = []);
        hi  = idsInto(voc[nm], tokens(row[cols[nm]], defs[nm].sep), out);
        src = out; lo = 0;
      }
      SC_SRC[i] = src; SC_LO[i] = lo; SC_HI[i] = hi;

      if (hasSel) {
        var hit = false, strs = voc[nm].strs;
        for (var k = lo; k < hi; k++) {
          if (chosen.has(strs[src[k]])) { hit = true; break; }
        }
        if (!hit) {
          fails++; failed = i;
          if (fails > 1) break;                           // counts for nobody; stop reading
        }
      }
    }

    if (fails > 1) return false;
    if (fails === 1) {                                    // AND across columns
      bump(cnt[names[failed]], SC_SRC[failed], SC_LO[failed], SC_HI[failed]);
      return false;
    }
    for (i = 0; i < names.length; i++) {
      bump(cnt[names[i]], SC_SRC[i], SC_LO[i], SC_HI[i]);
    }
    return true;
  }

  $.fn.dataTable.ext.search.push(facetPredicate);

  function scheduleDraw(w) {
    clearTimeout(PENDING[w.key]);
    PENDING[w.key] = setTimeout(function () {
      // resetPaging default: a filter change should return to the top, which is also what
      // Scroller does with the scroll position.
      w.api.draw();
    }, DRAW_DEBOUNCE_MS);
  }

  // ---- the open menu -------------------------------------------------------
  var EDGE = 8;    // keep this much clear of every viewport edge

  function closeOpen() {
    if (!OPEN) return;
    var w = OPEN;
    w.$menu.prop("hidden", true);
    w.$btn.attr("aria-expanded", "false");
    OPEN = null;
  }

  // The menu is fixed, so its coordinates are the button's viewport coordinates. Aligned
  // to the button's left edge, flipped to its right edge (and then clamped) rather than
  // being allowed off screen — the faceted columns sit far enough right in the SV table
  // for that to matter — and flipped above the header when there is no room below.
  function placeMenu(w) {
    var b    = w.$btn[0].getBoundingClientRect();
    var menu = w.$menu[0];
    var mw   = menu.offsetWidth, mh = menu.offsetHeight;

    var left = b.left;
    if (left + mw > window.innerWidth - EDGE) left = b.right - mw;
    left = Math.max(EDGE, Math.min(left, window.innerWidth - EDGE - mw));

    var top = b.bottom + 4;
    if (top + mh > window.innerHeight - EDGE) {
      var above = b.top - 4 - mh;
      top = above >= EDGE ? above
                          : Math.max(EDGE, window.innerHeight - EDGE - mh);
    }

    menu.style.left = Math.round(left) + "px";
    menu.style.top  = Math.round(top) + "px";
  }

  function openMenu(w) {
    if (OPEN === w) { closeOpen(); return; }
    closeOpen();
    w.$menu.prop("hidden", false);           // measurable only once it is not display:none
    w.$btn.attr("aria-expanded", "true");
    OPEN = w;
    placeMenu(w);
    // preventScroll: the menu is already in the viewport, and letting the browser scroll an
    // ancestor to "reveal" a focused field is how the table would jump under the reader.
    if (w.$find.length) w.$find[0].focus({ preventScroll: true });
  }

  // A fixed menu does not travel with the button, so it has to be told to. Scroll events do
  // not bubble but they do capture, so one document-level capturing listener catches the
  // page, the scroll body and the scroll head alike. Only ever one menu is open, and the
  // handler leaves immediately when there is none.
  function followButton() {
    if (!OPEN) return;
    var w = OPEN;
    var b = w.$btn[0].getBoundingClientRect();
    var clip = w.clip ? w.clip.getBoundingClientRect() : null;
    // Scrolled out from under its own column, or off the page entirely: a menu still
    // pointing at a column that is no longer there is worse than no menu.
    if ((b.width === 0 && b.height === 0) ||
        b.bottom < 0 || b.top > window.innerHeight ||
        (clip && (b.right <= clip.left || b.left >= clip.right))) {
      closeOpen();
      return;
    }
    placeMenu(w);
  }

  // ---- widget -------------------------------------------------------------
  // Two nodes, in two places: the button replaces the text box in the header cell, the menu
  // is appended to <body>. See note 1 at the top for why they cannot share a parent.
  function buttonMarkup(name) {
    return '<div class="facet">' +
             '<button type="button" class="facet__btn" aria-expanded="false"' +
             ' aria-haspopup="true" title="Filter ' + esc(name) + ' by value">' +
               '<span class="facet__label">All</span>' +
               '<span class="facet__caret" aria-hidden="true">&#9662;</span>' +
             '</button>' +
           '</div>';
  }

  function menuMarkup(name, def) {
    var withFind = def.values.length > SEARCH_THRESHOLD;
    var opts = def.values.map(function (v) {
      var tok = v[0], n = v[1];
      var isNone = (tok === null);
      var label = isNone ? "(none)" : tok;
      return '<label class="facet__opt' + (isNone ? " facet__opt--none" : "") + '"' +
             ' title="' + esc(label) + '">' +
             '<input type="checkbox"' +
             (isNone ? ' data-none="1"' : ' data-tok="' + esc(tok) + '"') + '>' +
             '<span class="facet__v">' + esc(label) + '</span>' +
             '<span class="facet__n">' + Number(n).toLocaleString() + '</span>' +
             '</label>';
    }).join("");

    return '<div class="facet__menu" role="group" aria-label="' + esc(name) +
           ' values" hidden>' +
             '<div class="facet__tools">' +
               (withFind ? '<input type="text" class="facet__find"' +
                           ' placeholder="Find value…" autocomplete="off">' : "") +
               '<button type="button" class="facet__all" title="Tick all shown">All' +
               '</button>' +
               '<button type="button" class="facet__none" title="Untick all shown">None' +
               '</button>' +
             '</div>' +
             '<div class="facet__list">' + opts + '</div>' +
             '<p class="facet__foot">counts follow the other active filters</p>' +
           '</div>';
  }

  function tokenOf($cb) {
    return $cb.data("none") ? null : String($cb.attr("data-tok"));
  }

  // A column with nothing ticked is deliberately absent from STATE rather than holding an
  // empty Set: the predicate reads STATE per row, 167k times per draw, and an absent key
  // is one property lookup where a Set would be a construction and a size test.
  var NO_SELECTION = new Set();

  function readSel(w) {
    return STATE[w.key][w.name] || NO_SELECTION;
  }

  function ensureSel(w) {
    var sel = STATE[w.key][w.name];
    if (!sel) { sel = new Set(); STATE[w.key][w.name] = sel; }
    return sel;
  }

  function pruneSel(w) {
    var sel = STATE[w.key][w.name];
    if (sel && sel.size === 0) delete STATE[w.key][w.name];
  }

  function trunc(s) {
    return s.length > LABEL_MAX ? s.slice(0, LABEL_MAX - 1) + "…" : s;
  }

  function refreshLabel(w) {
    var sel = readSel(w), n = sel.size, text, title = "";
    if (n === 0) {
      text = "All";
    } else if (n === 1) {
      var only = sel.values().next().value;
      var lab  = (only === null) ? "(none)" : only;
      text  = trunc(lab);
      title = lab;
    } else {
      text  = n + " selected";
      title = Array.from(sel).map(function (v) { return v === null ? "(none)" : v; })
                             .join(", ");
    }
    w.$label.text(text);
    w.$btn.attr("title", title || ("Filter " + w.name + " by value"));
    w.$container.toggleClass("facet--active", n > 0);
  }

  // Paint the counts accumulated by the last filtering pass. Raw DOM, and only where the
  // number actually changed: this runs on every draw, and the option rows live in the live
  // header DT re-measures.
  function flushCounts(key) {
    var cnt = COUNTS[key];
    if (!cnt) return;
    for (var i = 0; i < WIDGETS.length; i++) {
      var w = WIDGETS[i];
      if (w.key !== key) continue;
      var arr = cnt[w.name];
      if (!arr) continue;
      for (var j = 0; j < w.optIds.length; j++) {
        var n = arr[w.optIds[j]] || 0;
        if (w.optLast[j] === n) continue;
        w.optLast[j] = n;
        w.numNodes[j].textContent = n.toLocaleString();
        // The baseline R emitted is what the visible number no longer shows, so it moves
        // into the tooltip rather than being lost.
        w.optNodes[j].title = w.optLabel[j] + " — " + n.toLocaleString() + " of " +
                              w.optBase[j].toLocaleString() + " rows";
        w.optNodes[j].classList.toggle("facet__opt--empty", n === 0);
      }
    }
  }

  function anyActive() {
    return TABLES.some(function (t) {
      return Object.keys(STATE[t.key]).some(function (c) {
        return STATE[t.key][c] && STATE[t.key][c].size > 0;
      });
    });
  }

  // The report's own control bar owns the reset, because with scrollY = "400px" a ticked
  // filter can be scrolled out of sight and "the table is empty and I don't know why" is
  // the predictable question. Absent in an older template, hence the guard.
  function syncClearButton() {
    var el = document.getElementById("facet-clear");
    if (el) el.hidden = !anyActive();
  }

  function clearAll() {
    WIDGETS.forEach(function (w) {
      delete STATE[w.key][w.name];
      w.$menu.find('input[type="checkbox"]').prop("checked", false);
      refreshLabel(w);
    });
    syncClearButton();
    TABLES.forEach(function (t) {
      var node = t.elem();
      if (node && $.fn.dataTable.isDataTable(node)) $(node).DataTable().draw();
    });
  }

  function build(t, api, $td, name, def) {
    var $wrap  = $td.children("div").first();      // DT's div.form-group.has-feedback
    var $input = $wrap.children("input");

    // Hidden, not removed: DT holds a reference to this input, and the scroll-head
    // un-clipping handler is bound to $wrap itself. We never write to the input, so this
    // column's DataTables search stays "" and cannot collide with the facet filter.
    $input.hide().attr("tabindex", "-1").attr("aria-hidden", "true");
    $wrap.children("span.glyphicon").hide();

    var $container = $(buttonMarkup(name));
    $wrap.append($container);                      // inside $wrap, see note 1 at the top
    var $menu = $(menuMarkup(name, def)).appendTo(document.body);

    // What the button is clipped by, so an open menu can tell that its column has been
    // scrolled away. .dataTables_scrollHead exists only under scrollX/scrollY; the table
    // container is the honest fallback.
    var $cont = $(api.table().container());
    var $clip = $cont.find(".dataTables_scrollHead").first();

    var w = {
      key: t.key, name: name, api: api,
      clip: ($clip.length ? $clip : $cont)[0],
      $container: $container,
      $btn:   $container.find(".facet__btn"),
      $label: $container.find(".facet__label"),
      $menu:  $menu,
      $find:  $menu.find(".facet__find"),
      $list:  $menu.find(".facet__list")
    };

    // Option rows, in payload order, paired with their vocabulary id so the live count is
    // an array read. Cached as raw nodes: flushCounts() touches them on every draw.
    var voc = VOCAB[t.key][name];
    w.optNodes = w.$list.children(".facet__opt").toArray();
    w.numNodes = w.optNodes.map(function (el) { return el.querySelector(".facet__n"); });
    w.optIds   = def.values.map(function (v) { return idOf(voc, v[0]); });
    w.optBase  = def.values.map(function (v) { return Number(v[1]) || 0; });
    w.optLabel = def.values.map(function (v) { return v[0] === null ? "(none)" : v[0]; });
    w.optLast  = def.values.map(function () { return -1; });   // -1: nothing painted yet

    WIDGETS.push(w);

    w.$btn.on("click", function (e) {
      e.stopPropagation();
      openMenu(w);
    });

    // Delegated within our own node — never from `document`, which the per-draw header
    // clone would also match.
    w.$menu.on("change", 'input[type="checkbox"]', function () {
      var sel = ensureSel(w), tok = tokenOf($(this));
      if (this.checked) sel.add(tok); else sel.delete(tok);
      pruneSel(w);
      refreshLabel(w);
      syncClearButton();
      scheduleDraw(w);
    });

    w.$menu.find(".facet__all, .facet__none").on("click", function () {
      var on  = $(this).hasClass("facet__all");
      var sel = ensureSel(w);
      // Only what is currently shown, which is what the button titles say — once the find
      // box has narrowed the list, acting on the hidden values too would be a surprise.
      w.$list.children(".facet__opt").not(".facet__opt--hidden")
        .find('input[type="checkbox"]').each(function () {
          this.checked = on;
          var tok = tokenOf($(this));
          if (on) sel.add(tok); else sel.delete(tok);
        });
      pruneSel(w);
      refreshLabel(w);
      syncClearButton();
      scheduleDraw(w);
    });

    if (w.$find.length) {
      w.$find.on("input", function () {
        var q = this.value.trim().toLowerCase();
        w.$list.children(".facet__opt").each(function () {
          var txt = $(this).find(".facet__v").text().toLowerCase();
          $(this).toggleClass("facet__opt--hidden", q.length > 0 && txt.indexOf(q) < 0);
        });
      });
    }

    // Bound to both nodes: they are no longer in the same subtree, and Escape has to work
    // from the find box as much as from the button.
    w.$container.add(w.$menu).on("keydown", function (e) {
      if (e.key === "Escape" || e.keyCode === 27) {
        closeOpen();
        w.$btn.trigger("focus");
      }
    });

    refreshLabel(w);
    return w;
  }

  // ---- attach -------------------------------------------------------------
  // Note 5 at the top: whatever order the report's scripts happen to be emitted in, the
  // counts only mean "what is left after everything else" if we filter last.
  function moveToEnd() {
    var list = $.fn.dataTable.ext.search;
    var at = list.indexOf(facetPredicate);
    if (at >= 0 && at !== list.length - 1) {
      list.splice(at, 1);
      list.push(facetPredicate);
    }
  }

  function attach(t) {
    if (ATTACHED[t.key]) return;
    var node = t.elem();
    if (!node || !$.fn.dataTable.isDataTable(node)) return;
    var defs = t.facets(), cols = t.cols();
    if (!defs || !cols) return;

    var api    = $(node).DataTable();
    var $head  = $(api.table().header());          // the live thead, wherever DT moved it
    var $rows  = $head.children("tr");
    if ($rows.length < 2) return;                  // no filter row: nothing to replace
    ATTACHED[t.key] = true;
    ensureTable(t);
    moveToEnd();

    var $label  = $rows.first().children("th,td");
    var $filter = $rows.last().children("td");
    var visIdx  = api.columns(":visible").indexes().toArray();

    if ($filter.length !== visIdx.length) {
      console.warn("facet_filter: " + t.key + " filter row has " + $filter.length +
                   " cells for " + visIdx.length +
                   " visible columns — leaving the text filters alone.");
      return;
    }

    // R drops a column it cannot offer a dropdown for — fewer than two distinct values, or
    // gone from the frame entirely — and knitr swallows its message(), so this is where that
    // decision becomes visible. console.info rather than warn: for `callers` on the VEP text
    // path and the SV `caller` column it is the expected outcome, and warnings here are
    // reserved for something actually being wrong.
    var requested = t.requested() || [];
    var missing = requested.filter(function (n) { return !(n in defs); });
    if (missing.length) {
      console.info("facet_filter: no value list for " + t.key + " column(s) " +
                   missing.join(", ") + " (fewer than " +
                   "2 distinct values, or not in the table) — they keep their text filter.");
    }

    Object.keys(defs).forEach(function (name) {
      var idx = cols[name];
      if (idx === undefined) {
        console.warn("facet_filter: '" + name + "' is not in " + t.key.toUpperCase() +
                     "_COLS — was the column renamed?");
        return;
      }
      var pos = visIdx.indexOf(idx);
      if (pos < 0) {
        console.warn("facet_filter: '" + name + "' is a hidden column.");
        return;
      }
      // Decisive cross-check on the translation: the label above the cell we picked must be
      // the column we think it is. Attaching an svtype dropdown to the locus column would
      // otherwise look like a filter that works and filters the wrong thing.
      var got = $.trim($label.eq(pos).text());
      if (got !== name) {
        console.warn("facet_filter: header mismatch for '" + name + "' at visible column " +
                     pos + " (found '" + got + "') — skipping.");
        return;
      }
      var $td = $filter.eq(pos);
      if ($td.attr("data-type") !== "character") {
        // A factor/logical column gets DT's own selectize control; two widgets in one cell
        // would fight over the same input.
        console.warn("facet_filter: '" + name + "' has data-type '" +
                     $td.attr("data-type") + "' — leaving DT's own filter in place.");
        return;
      }
      build(t, api, $td, name, defs[name]);
    });

    syncClearButton();

    // Speed only — the predicate splits live until this exists, and stays correct either
    // way. Guarded because a failure here must cost performance, not the filter.
    try {
      buildIndex(t, api);
    } catch (err) {
      console.warn("facet_filter: could not index " + t.key +
                   " tokens, falling back to splitting per draw.", err);
    }

    // DataTables fires `search` at the end of _fnFilterComplete, i.e. immediately after the
    // pass that filled the counts — no assumption about draw ordering needed. A pass that
    // fires it without ever calling the predicate filtered everything out upstream, and the
    // counts are then genuinely zero.
    $(node).on("search.dt", function (e) {
      if (e.target !== node) return;
      if (!TOUCHED[t.key]) resetCounts(t.key);
      TOUCHED[t.key] = false;
      flushCounts(t.key);
    });

    // The numbers rendered by R are the unfiltered baseline, and --gene-panel means the
    // table can already be filtered on load. One forced pass now, so what is on screen is
    // what the predicate says rather than what the payload said.
    api.draw(false);
  }

  // DT wires its own filter-row handlers after $().DataTable() returns, i.e. after init.dt
  // has fired, and one of those is the handler that un-clips the scroll head for an open
  // menu. Defer a tick so it exists before we touch the cell. The load sweep is the
  // idempotent fallback for a table that initialised before this script was parsed.
  $(document).on("init.dt", function (e, settings) {
    var t = tableFor(settings.nTable);
    if (t) setTimeout(function () { attach(t); }, 0);
  });

  $(window).on("load", function () {
    TABLES.forEach(function (t) { if (t.elem()) attach(t); });
    var clear = document.getElementById("facet-clear");
    if (clear) clear.addEventListener("click", clearAll);
    syncClearButton();
  });

  // The one outside-click handler. Bound to document on purpose — it tests containment
  // against the open widget's own two nodes, not against a widget class, so the header's
  // sizing clone is irrelevant to it. Both nodes, because the menu is in <body> and the
  // button in the header: a click inside the menu is not inside the button's container.
  $(document).on("mousedown", function (e) {
    if (!OPEN) return;
    if (OPEN.$container[0].contains(e.target) || OPEN.$menu[0].contains(e.target)) return;
    closeOpen();
  });

  // A fixed menu has to be re-placed whenever anything under it moves. Capturing, so the
  // scroll body's and the scroll head's own scroll events reach it.
  document.addEventListener("scroll", followButton, true);
  window.addEventListener("resize", followButton);

  // Which values are ticked, for inspection from the browser console and for
  // tests/js/test-facet-predicate.js — the ticked state is otherwise unreachable without a
  // DOM. Same reason window.svPanelState is exposed in per_sample.qmd: one shared object
  // beats re-deriving the state from the rendered header.
  window.facetFilterState = STATE;

  // The live counts of the last filtering pass, as { table: { column: [[token, n], ...] } }
  // in payload order, for the same two audiences. A function, because the internal form is
  // integer-keyed and only meaningful beside its vocabulary.
  window.facetFilterCounts = function () {
    var out = {};
    Object.keys(COUNTS).forEach(function (key) {
      var cnt = COUNTS[key], voc = VOCAB[key];
      out[key] = {};
      Object.keys(cnt).forEach(function (nm) {
        // Over the vocabulary, not over the count array: a value that matched no row is a
        // 0, not a gap, and the two lengths differ whenever a token was added after the
        // last reset.
        var strs = voc[nm].strs, pairs = [];
        for (var id = 0; id < strs.length; id++) pairs.push([strs[id], cnt[nm][id] || 0]);
        out[key][nm] = pairs;
      });
    });
    return out;
  };
})();
