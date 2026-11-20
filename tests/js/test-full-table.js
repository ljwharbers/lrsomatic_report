// Decode / download / unpack tests for the embedded full small-variant table
// (assets/js/full_table.js), on a sidecar written by write_full_variant_table()
// (tests/fixtures/full_table/fixture_small_variants.tsv.gz: quoted fields, NA, logicals, "&"/">").
//
// Plain `node`, no dependencies. Needs Node >= 18 for Blob/DecompressionStream/Response.
//
//   node tests/js/test-full-table.js

const fs = require("fs");
const path = require("path");
const zlib = require("zlib");
const ft = require(path.join(__dirname, "..", "..", "assets", "js", "full_table.js"));

let failed = 0;
const pending = [];
function check(name, got, want) {
  const g = JSON.stringify(got), w = JSON.stringify(want);
  if (g !== w) {
    failed++;
    console.error("FAIL " + name + ": got " + g + ", want " + w);
  }
}

const gz = fs.readFileSync(path.join(__dirname, "..", "fixtures", "full_table",
                                     "fixture_small_variants.tsv.gz"));
const tsv = zlib.gunzipSync(gz).toString("utf8");

// base64 -> bytes is byte-exact, and tolerates the line breaks/whitespace an encoder may add
const b64 = gz.toString("base64");
check("decode", Buffer.from(ft.b64ToBytes(b64)).equals(gz), true);
check("decode, wrapped", Buffer.from(ft.b64ToBytes(b64.replace(/(.{60})/g, "$1\n  "))).equals(gz), true);

// Download hands back the same bytes under the sidecar's name, then frees the URL
(function () {
  const made = [], clicked = [], revoked = [];
  const doc = {
    body: { appendChild: (a) => made.push(a), removeChild: () => {} },
    createElement: () => ({ style: {}, click() { clicked.push(this.download); } })
  };
  let blob = null;
  const urlApi = { createObjectURL: (b) => { blob = b; return "blob:x"; },
                   revokeObjectURL: (u) => revoked.push(u) };
  const realSetTimeout = global.setTimeout;
  global.setTimeout = (fn) => fn();
  try {
    ft.downloadBytes(ft.b64ToBytes(b64), "S1_small_variants.tsv.gz", doc, urlApi);
  } finally {
    global.setTimeout = realSetTimeout;
  }
  check("download filename", clicked, ["S1_small_variants.tsv.gz"]);
  check("download href", made[0].href, "blob:x");
  check("download type", blob.type, "application/gzip");
  check("download revoked", revoked, ["blob:x"]);
  pending.push(blob.arrayBuffer().then((ab) => check("download bytes", Buffer.from(ab).equals(gz), true)));
})();

// Feature detection: no DecompressionStream, no "show all"
check("canDecompress missing", ft.canDecompress({ Blob: function () {}, Response: function () {} }), false);
check("canDecompress present", ft.canDecompress({ DecompressionStream: function () {},
                                                  Blob: function () {}, Response: function () {} }), true);

// fwrite's quote = "auto": doubled quotes, tab and newline inside quotes, "" for an empty string
const parsed = ft.parseTsv(tsv);
check("header", parsed.header,
      ["symbol", "chrom", "pos", "ref", "alt", "consequence", "impact", "hgvsc", "vaf", "dp",
       "note", "flag"]);
check("row count", parsed.rows.length, 4);
check("quoted quote", parsed.rows[1][10], 'say "hi"');
check("quoted tab", parsed.rows[2][10], "tab\there");
check("quoted newline", parsed.rows[3][10], "two\nlines");
check("empty string", parsed.rows[3][5], "");
check("unquoted fast path", ft.parseTsv("a\tb\r\n1\t\n").rows, [["1", ""]]);

// Rows land by SNV_COLS name; VAF% is vaf x 100 at 1 dp; NA/"" -> null; strings escaped as DT does
const spec = [["symbol", "symbol", "string", 1, null], ["pos", "pos", "number", 1, null],
              ["VAF%", "vaf", "number", 100, 1], ["dp", "dp", "number", 1, null],
              ["hgvsc", "hgvsc", "string", 1, null], ["flag", "flag", "logical", 1, null],
              ["consequence", "consequence", "string", 1, null],
              ["absent", "absent", "string", 1, null]];
const cols = { symbol: 0, pos: 1, "VAF%": 2, dp: 3, hgvsc: 4, flag: 5, consequence: 6, absent: 7 };
const rows = ft.toTableRows(parsed, spec, cols);
check("row 0", rows[0], ["TP53", 7675088, 45.7, 40, "c.524G&gt;A", true,
                         "missense_variant", null]);
check("row 1", rows[1], ["KRAS", 25245350, 12.5, null, "c.35G&gt;T", false,
                         "missense_variant,splice_region_variant", null]);
check("row 2", rows[2], [null, 100, null, 12, null, null, "intron_variant", null]);
check("row 3", rows[3], ["A&amp;B", 5, 0, 3, "c.1_2insT", true, null, null]);
check("range", ft.numericRange(rows, 1), [5, 25245350]);
check("range, all null", ft.numericRange(rows, 7), null);

// VAF% must round exactly as R's round(vaf * 100, 1) did for the embedded rows; expected values
// from R 4.4 (31.25 ties to even, 21.05000000000000071 is an R tie, 8.75 rounds up to even)
const vafs = [0.3125, 0.2105, 0.0875, 0.5455, 0.0625, 0.4567, 1e-4, 0.999999];
check("roundR matches R", vafs.map((v) => ft.roundR(v * 100, 1)),
      [31.2, 21, 8.8, 54.5, 6.2, 45.7, 0, 100]);
check("roundR negative", ft.roundR(-31.25, 1), -31.2);

// The browser path: gunzip the decoded bytes in-process
if (ft.canDecompress(globalThis)) {
  pending.push(ft.gunzipToText(ft.b64ToBytes(b64)).then((text) => check("gunzip", text, tsv)));
} else {
  console.log("skip: this node has no DecompressionStream");
}

Promise.all(pending).catch((e) => { failed++; console.error("FAIL " + e); }).then(() => {
  if (failed) { console.error(failed + " check(s) failed"); process.exit(1); }
  console.log("test-full-table: all checks passed (" + pending.length + " async)");
});
