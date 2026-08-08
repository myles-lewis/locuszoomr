/* Client-side gene track re-layout for locuszoomr.
 *
 * Injected into plotly widgets by add_genetrack_relayout(). Re-packs the gene
 * annotation panel against the visible window whenever the user zooms or pans.
 * Plain ES5 on purpose: this file ships as package source with no build step.
 */
(function (global) {
  'use strict';

  var LZR = global.LZR || {};

  /* Greedy first-fit packing. items: [{min, max}] in priority order.
   * Returns an array of 1-based row numbers, mirroring R's pack_rows(). */
  LZR.packRows = function (items) {
    var rows = [];
    var out = new Array(items.length);
    for (var i = 0; i < items.length; i++) {
      var it = items[i];
      var placed = false;
      for (var j = 0; j < rows.length && !placed; j++) {
        var clash = false;
        for (var k = 0; k < rows[j].length; k++) {
          var h = rows[j][k];
          if (it.min < h.max && it.max > h.min) { clash = true; break; }
        }
        if (!clash) { rows[j].push(it); out[i] = j + 1; placed = true; }
      }
      if (!placed) { rows.push([it]); out[i] = rows.length; }
    }
    return out;
  };

  /* Convert a trace axis id ("x", "x2") to its layout key ("xaxis", "xaxis2"). */
  LZR.axisKey = function (id) {
    return id.charAt(0) + 'axis' + id.slice(1);
  };

  var measureCtx = null;
  LZR.measure = function (text, fontPx) {
    if (measureCtx === null) {
      measureCtx = document.createElement('canvas').getContext('2d');
    }
    measureCtx.font = fontPx + 'px "Open Sans", verdana, arial, sans-serif';
    return measureCtx.measureText(text).width;
  };

  /* Pack the genes visible in [x0, x1] into rows.
   * Returns {rows: [...], nrow: n, hidden: n} where rows entries carry the
   * gene index, its assigned row, and its label anchor. */
  LZR.layout = function (P, x0, x1, pxPerData) {
    var span = x1 - x0;
    var gap = span * P.cfg.gapFrac;
    var margin = span * 0.1;
    var lo = x0 - margin, hi = x1 + margin;

    var items = [];
    for (var i = 0; i < P.tx.length; i++) {
      var g = P.tx[i];
      if (g.end < lo || g.start > hi) continue;
      /* mirror R: strwidth(paste0("--", gene_name)) */
      var halfw = 0;
      if (g.label !== '') {
        halfw = LZR.measure('--' + g.name, P.cfg.fontSizePx) / pxPerData / 2;
      }
      var mid = (g.start + g.end) / 2;

      /* Edge-label clamping, mirroring mapRow() in R/genetracks.R. A gene
       * straddling the viewport edge would otherwise anchor its label at its
       * own midpoint, which can sit off-screen, so the bar renders but the
       * name does not. Pull the anchor just inside the boundary instead --
       * but only when the gene is wide enough to hold the whole label, else
       * the text would float beyond the gene it belongs to.
       *
       * The conditions and the sequential left-then-right order match R (a
       * gene wider than the view can satisfy both, and the right-hand test
       * sees the already-clamped anchor), but the boundary does NOT. R
       * clamps to xlim widened by 4%, which works on a graphics device
       * because text may render into the plot margin. Plotly clips hard at
       * the axis range, so an anchor placed at `x0 - 0.04 * span + halfw`
       * is still off-screen whenever halfw is smaller than that overshoot,
       * and only the tail of the label (the strand arrow) shows. Clamp to
       * the visible window itself so the whole label lands inside it.
       *
       * Unlabelled genes are skipped: halfw is 0 for them, which would
       * satisfy the conditions trivially and shift the packing footprint
       * for no visible gain. */
      if (halfw > 0) {
        var gwFull = halfw * 2;
        if ((mid - halfw) < x0 && (x0 + gwFull) < g.end) mid = x0 + halfw;
        if ((mid + halfw) > x1 && (x1 - gwFull) > g.start) mid = x1 - halfw;
      }
      items.push({
        i: i, mid: mid, priority: g.priority,
        min: Math.min(g.start, g.end, mid - halfw) - gap / 2,
        max: Math.max(g.start, g.end, mid + halfw) + gap / 2
      });
    }
    items.sort(function (a, b) { return a.priority - b.priority; });

    var assigned = LZR.packRows(items);
    var nrow = 0;
    for (var k = 0; k < assigned.length; k++) {
      items[k].row = assigned[k];
      if (assigned[k] > nrow) nrow = assigned[k];
    }

    var kept = [], hidden = 0;
    for (var m = 0; m < items.length; m++) {
      if (items[m].row <= P.cfg.maxrows) kept.push(items[m]); else hidden++;
    }
    return {rows: kept, nrow: nrow, hidden: hidden};
  };

  LZR.attach = function (el, P) {
    var gd = el;
    var dead = false;

    function fail(err) {
      dead = true;
      if (global.console) global.console.warn('locuszoomr re-layout disabled:', err);
    }

    /* Everything below can throw synchronously (e.g. a malformed payload),
     * and an uncaught throw here would escape the onRender hook entirely,
     * abort remaining htmlwidgets afterRender hooks on this element, and
     * leave no re-layout listener attached. Fail closed instead: warn once
     * and leave the widget exactly as R rendered it. */
    try {
      var Plotly = global.Plotly;
      var xkey = LZR.axisKey(P.idx.xaxis);
      var yref = P.idx.yaxis;

      /* Shapes that are not ours, captured once. Every update rebuilds the
       * array as base.concat(exons), so stored indices never go stale.
       * shapeIdx serialises as a bare number, not a length-1 array, when the
       * gene panel has exactly one shape (jsonlite auto_unbox); normalise to
       * an array before relying on .indexOf(). Use `== null` rather than
       * `||` here: 0 is a valid (and the most common) shape index, and 0 is
       * falsy, so `P.idx.shapeIdx || []` would silently discard exactly the
       * single-shape-at-index-0 case this normalisation exists for. */
      var idxs = [].concat(P.idx.shapeIdx == null ? [] : P.idx.shapeIdx);
      var all = (gd.layout.shapes || []);
      var base = [];
      for (var i = 0; i < all.length; i++) {
        if (idxs.indexOf(i) === -1) base.push(all[i]);
      }
      var baseAnn = (gd.layout.annotations || []).slice();

      var busy = false, timer = null, pending = false, lastKey = null;

      /* Plotly.relayout(gd, {shapes: ..., annotations: ...}) inside apply()
       * below re-emits 'plotly_relayout' on this same gd from inside its own
       * .then() chain, while `busy` is still true. That self-event now takes
       * the `pending` path (see the listener below) rather than being
       * dropped outright, so *something* has to stop it from re-triggering
       * apply() forever. This function is that something: it recognises an
       * event as self-emitted, structurally, by checking that every key on
       * the update object is one we ourselves write. Do not rely on the
       * no-op range/length guard inside apply() for this instead — that
       * guard exists for a different, unrelated reason (see the comment at
       * its call site) and is not a safe substitute: relaxing or removing it
       * must not be able to reopen this loop. */
      function isSelfUpdate(upd) {
        if (!upd) return false;
        var keys = Object.keys(upd);
        if (keys.length === 0) return false;
        for (var k = 0; k < keys.length; k++) {
          if (keys[k] !== 'shapes' && keys[k] !== 'annotations') return false;
        }
        return true;
      }

      function schedule() {
        if (timer) global.clearTimeout(timer);
        timer = global.setTimeout(function () {
          try { apply(); } catch (err) { fail(err); }
        }, 100);
      }

      function apply() {
        var ax = gd._fullLayout[xkey];
        var rng = ax.range;
        var len = ax._length;
        /* Hidden-container guard: inside a Bootstrap tabset, flexdashboard
         * page, or any display:none ancestor, Plotly.Plots.resize() (fired
         * on tab switch) can emit 'plotly_relayout' while this axis has
         * zero (or as-yet-unset) pixel length. Neither failure mode throws,
         * so fail()/catch() never sees it: len === 0 makes pxPerData 0 and
         * every gene's footprint [-Infinity, +Infinity], so everything
         * clashes and collapses to one gene per row (silently dropping the
         * rest past maxrows); len === undefined makes pxPerData NaN, so no
         * clash comparison is ever true and every gene lands on row 1,
         * stacked on top of each other. Bail out before any of that,
         * *before* `busy` is set below, so neither it nor `pending` is left
         * stuck by this early return. */
        if (!(len > 0) || !isFinite(rng[0]) || !isFinite(rng[1])) return;
        /* No-op guard: dragmode toggles and legend clicks fire
         * plotly_relayout without changing the x range or the axis's
         * rendered pixel length, so skip the re-pack in that case. A window
         * resize is NOT such a no-op: it changes `ax._length` while leaving
         * `rng` unchanged, and pxPerData (hence every label-width
         * calculation below) is derived from `_length`, so the cache key
         * must include it too, not just the range. This guard is purely a
         * performance/cosmetic optimisation — it is NOT what stops the
         * self-triggered update loop described above; that termination is
         * structural, via isSelfUpdate() on the listener, and must keep
         * working even if this guard is later relaxed or removed. */
        if (lastKey !== null &&
            rng[0] === lastKey[0] && rng[1] === lastKey[1] && len === lastKey[2]) {
          return;
        }
        lastKey = [rng[0], rng[1], len];
        var pxPerData = len / (rng[1] - rng[0]);
        var res = LZR.layout(P, rng[0], rng[1], pxPerData);

        var lx = [], ly = [], lt = [], tx = [], ty = [], tt = [];
        var shapes = base.slice();
        for (var i2 = 0; i2 < res.rows.length; i2++) {
          var it = res.rows[i2], g = P.tx[it.i], y = -it.row;
          lx.push(g.start, g.end, null);
          ly.push(y, y, null);
          lt.push(g.hover, g.hover, null);
          if (g.label !== '') {
            tx.push(it.mid); ty.push(y + 0.35); tt.push(g.label);
          }
        }

        if (P.cfg.showExons === false) {
          /* Mirror R/genetrack_ly.R:178-185: one rect per gene, not per exon. */
          for (var g2 = 0; g2 < res.rows.length; g2++) {
            var itg = res.rows[g2], gg = P.tx[itg.i], rr = itg.row;
            shapes.push({
              type: 'rect', fillcolor: P.cfg.geneCol,
              line: {color: P.cfg.exonBorder, width: 1},
              x0: gg.start, x1: gg.end, xref: P.idx.xaxis,
              y0: -rr - 0.15, y1: -rr + 0.15, yref: yref
            });
          }
        } else {
          for (var j = 0; j < P.ex.length; j++) {
            var e = P.ex[j];
            var r = null;
            for (var k = 0; k < res.rows.length; k++) {
              if (res.rows[k].i === e.gene_idx) { r = res.rows[k].row; break; }
            }
            if (r === null) continue;
            shapes.push({
              type: 'rect', fillcolor: P.cfg.exonCol,
              line: {color: P.cfg.exonBorder, width: 0.5},
              x0: e.start, x1: e.end, xref: P.idx.xaxis,
              y0: -r - 0.15, y1: -r + 0.15, yref: yref
            });
          }
        }

        /* KNOWN UNFIXED LIMITATION: this annotation branch only executes
         * inside apply(), which itself only runs in response to a
         * 'plotly_relayout' event (see gd.on(...) below and the initial
         * onRender wiring in R/genetrack_relayout.R). It never runs on the
         * INITIAL render — the widget is handed to the browser already
         * built by R, with no relayout event to trigger a re-pack. That is
         * exactly the view where truncation (rows beyond maxrows) is most
         * likely, since R packed at a possibly-too-narrow `width`. So on
         * first paint there is no on-plot "N genes not shown" annotation
         * at all; the only truncation signal available to the user at that
         * point is the console message() emitted R-side in
         * R/genetrack_ly.R ("N tracks needed to show all genes"), which is
         * invisible unless the browser/R console is open. */
        var ann = baseAnn.slice();
        if (res.rows.length === 0) {
          ann.push(LZR.note('No genes in view', P));
        } else if (res.hidden > 0) {
          ann.push(LZR.note(res.hidden + ' genes not shown — zoom in', P));
        }

        busy = true;
        Plotly.restyle(gd, {x: [lx], y: [ly], text: [lt]}, [P.idx.lineTrace])
          .then(function () {
            return Plotly.restyle(gd, {x: [tx], y: [ty], text: [tt]},
                                  [P.idx.labelTrace]);
          })
          .then(function () {
            return Plotly.relayout(gd, {shapes: shapes, annotations: ann});
          })
          .then(function () {
            busy = false;
            /* A *genuine* relayout (not our own Plotly.relayout(shapes/
             * annotations) call above, which the listener below already
             * filters out via isSelfUpdate before pending is ever set) can
             * arrive while this apply() was in flight; the three-restyle
             * chain takes long enough (150-300ms on a dense locus) that a
             * dragmode="pan" user can easily pan again before it settles.
             * Re-run once more against the now-current range instead of
             * leaving the panel packed for the stale window. This does NOT
             * on its own risk an infinite loop: Plotly.relayout's own
             * self-emitted event is excluded upstream by isSelfUpdate(), so
             * `pending` only ever becomes true for a real user-driven event. */
            if (pending) { pending = false; schedule(); }
          })
          .catch(function (err) { busy = false; fail(err); });
      }

      gd.on('plotly_relayout', function (upd) {
        if (dead) return;
        /* Structural loop-breaker for the self-triggering update described
         * above isSelfUpdate()'s definition: without this, our own
         * Plotly.relayout({shapes, annotations}) call would re-enter here
         * via `pending`/schedule() and run forever. */
        if (isSelfUpdate(upd)) return;
        if (busy) { pending = true; return; }
        schedule();
      });
    } catch (err) {
      fail(err);
    }
  };

  LZR.note = function (text, P) {
    return {
      text: text, showarrow: false,
      xref: P.idx.xaxis + ' domain', yref: P.idx.yaxis + ' domain',
      x: 1, y: 0, xanchor: 'right', yanchor: 'bottom',
      font: {size: P.cfg.fontSizePx * 0.9, color: '#888888'}
    };
  };

  global.LZR = LZR;
})(typeof window !== 'undefined' ? window : global);
