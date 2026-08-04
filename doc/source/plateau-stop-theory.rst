.. _plateauStopTheory:

Plateau-based MIP stopping: theoretical description
*****************************************************

This page describes the algorithm behind the :ref:`plateauStop` /
:ref:`plateauStopBound` features (options ``mip:plateautime``,
``mip:plateauabstol``, ``mip:plateaureltol``, ``mip:plateauwarmup``,
``mip:plateauwarmuprelgap``, ``mip:plateauwarmupabsgap``,
``mip:plateauabsgaptol``, ``mip:plateaurelgaptol``), implemented once
in `MIPBackend` (``include/mp/backend-mip.h``).
For the practical, option-by-option description see :ref:`plateauStop`; 
for the driver-implementation steps see :ref:`implement-standard-features`.


Motivation
==========

A MIP solver's branch-and-bound search produces two monotonic
sequences over time:

* the **incumbent objective**, which only ever improves (a new
  incumbent is by definition better than the previous one), and
* the **best (dual) bound**, which only ever tightens towards the
  incumbent, and whose distance from the incumbent is the **MIP gap**.

Classic stopping criteria act on the *values* of these sequences: a
time limit bounds wall-clock time regardless of progress, and a gap
limit stops once the two sequences are close enough to each other.
Neither criterion looks at the *shape* of the sequences over time.

In practice, a search often reaches a point where further node
exploration keeps technically improving the incumbent and/or narrowing
the gap, but by amounts so small, so infrequently, that continuing is
not worth the wall-clock cost. This is the *plateau*: the search
hasn't necessarily converged (the gap may still be large), but it has
stopped making an amount of progress a user considers worthwhile per
unit of time. The feature described here detects that condition and
proactively terminates the solve, independently of and in addition to
time/gap limits.


Formal model: the incumbent channel
=====================================

The incumbent channel is always active once ``mip:plateautime`` is set
(:ref:`plateauStop`, no extra solver support needed beyond the base
feature). The driver reports a sequence of incumbent objective values
:math:`v_1, v_2, \dots` at times :math:`t_1 < t_2 < \dots`, each
:math:`v_k` strictly better than the last -- guaranteed automatically,
because a native "new incumbent" callback only fires on genuine
improvements. The algorithm therefore never needs to know the
optimization sense (minimize/maximize) and works purely with absolute
differences.

The state consists of a reference (baseline) value :math:`\hat v`, the
time of the last *qualifying* update :math:`t_{\text{last}}`, and the
time the feature was armed (solve start) :math:`t_0`.

**Progress predicate.** A report :math:`v_k` counts as progress
against the current baseline :math:`\hat v` iff

.. math::

  P(\hat v, v_k) \;=\;
    \big|\hat v - v_k\big| > \varepsilon_{\text{abs}}
    \;\lor\;
    \frac{\big|\hat v - v_k\big|}{\big|\hat v\big|} > \varepsilon_{\text{rel}}

with the usual convention at :math:`\hat v = 0`: the relative term is
:math:`+\infty` if :math:`v_k \ne \hat v` (i.e. any nonzero change from
a zero baseline counts as progress under the relative test), and
:math:`0` otherwise. :math:`\varepsilon_{\text{abs}}` and
:math:`\varepsilon_{\text{rel}}` are ``mip:plateauabstol`` and
``mip:plateaureltol`` respectively, both defaulting to **0** -- which,
for this channel, means "any improvement, however small, counts",
matching the same "either absolute or relative" convention already
used by ``mip:gap``/``mip:gapabs``. (The gap channel below reuses this
same predicate shape, but with a different, opt-in, meaning for a zero
tolerance -- see below.)

**Update rule.** On each report :math:`v_k` at time :math:`t_k`:

.. math::

  \text{if } P(\hat v, v_k):\quad
    \hat v \leftarrow v_k, \quad t_{\text{last}} \leftarrow t_k

  \text{otherwise:}\quad \hat v, \, t_{\text{last}} \text{ unchanged}

Crucially, the baseline :math:`\hat v` only advances when a report
clears the tolerance. A report that does not clear it is *absorbed*:
it is simply discarded, without moving the reference point and
without resetting the clock. This matters because it makes the
predicate compare against the last point that actually counted as
progress, not against the immediately preceding report. Comparing
step-by-step instead (advancing :math:`\hat v` on every report
regardless of whether it qualified) would let an arbitrarily long
sequence of individually-negligible micro-improvements silently
accumulate real, tolerance-clearing progress while the clock keeps
being reset on every single one of them -- the search would never be
recognized as having plateaued even though, from the baseline's point
of view, nothing meaningful has happened.

**Stopping condition.** At any time :math:`t \ge t_0`, the search
should stop iff both:

.. math::

  \text{warmup done at } t
  \qquad\text{and}\qquad
  t - t_{\text{last}} \;\ge\; T_{\text{plateau}}

where :math:`T_{\text{plateau}}` is ``mip:plateautime`` and "warmup
done" is defined below. The condition is evaluated lazily, on every
report (or on any explicit check -- see :ref:`the periodic-callback
utility <plateauPeriodicUtility>` below), rather than on a separate
timer thread: since reports arrive on the solver's own callback
thread, the elapsed time is simply measured against a steady clock
read at each call, and termination (e.g. ``GRBterminate``) is
requested from within that same callback.


Two channels, asymmetric zero-tolerance semantics
====================================================

``mip:plateauabsgaptol`` / ``mip:plateaurelgaptol`` (gated behind
:ref:`plateauStopBound`, independent of each other) add a second,
optional *gap* channel, tracked in addition to -- not instead of --
the incumbent channel above:

.. math::

  P_{\text{absgap}}(\hat g_a, a_k) = \big|\hat g_a - a_k\big| > \varepsilon_{\text{absgap}}
  \qquad
  P_{\text{relgap}}(\hat g_r, r_k) = \big|\hat g_r - r_k\big| > \varepsilon_{\text{relgap}}

where :math:`a_k`/:math:`r_k` are the reported absolute/relative MIP
gap and :math:`\varepsilon_{\text{absgap}}`/:math:`\varepsilon_{\text{relgap}}`
are ``mip:plateauabsgaptol``/``mip:plateaurelgaptol``. Each sub-channel
uses the *same* update rule as the incumbent channel (baseline only
advances on a qualifying change) and shares the *same* clock
:math:`t_{\text{last}}` -- whichever channel (incumbent, absolute MIP gap,
relative MIP  gap) last produced a qualifying report resets it, so enabling more
channels can only make the plateau *harder* to reach, not easier: it
never introduces a new way to stop early, only new ways to keep the
timer alive.

The critical difference from the incumbent channel is what a zero
tolerance means. For ``mip:plateauabstol``/``mip:plateaureltol``, `0`
means "maximally sensitive: any change counts". 

For the gap sub-channels, `0` means **disabled**: that
sub-channel is not evaluated at all, and never contributes to
progress. 

With both gap sub-channels active alongside the incumbent channel, the
combined progress predicate is a plain logical OR:

.. math::

  \text{progress at } t_k \;=\;
    P(\hat v, v_k)
    \;\lor\;
    \big[\varepsilon_{\text{absgap}} > 0 \land P_{\text{absgap}}(\hat g_a, a_k)\big]
    \;\lor\;
    \big[\varepsilon_{\text{relgap}} > 0 \land P_{\text{relgap}}(\hat g_r, r_k)\big]



Role of the warm-up period
===========================

"Warmup done" (used in the stopping condition above) is true as soon
as *either* of two independent conditions holds:

.. math::

  t - t_0 \;\ge\; T_{\text{warmup}}
  \qquad\text{or}\qquad
  \big(\varepsilon_{\text{wgap,rel}} > 0 \land r \le \varepsilon_{\text{wgap,rel}}\big)
  \;\lor\;
  \big(\varepsilon_{\text{wgap,abs}} > 0 \land a \le \varepsilon_{\text{wgap,abs}}\big)

where :math:`T_{\text{warmup}}` is ``mip:plateauwarmup``,
:math:`\varepsilon_{\text{wgap,rel}}`/:math:`\varepsilon_{\text{wgap,abs}}`
are ``mip:plateauwarmuprelgap``/``mip:plateauwarmupabsgap`` (gated
behind :ref:`plateauStopBound`, like the gap-progress options
themselves -- an incumbent-only driver never has a real gap value to
offer the check below, so these two options are not even registered
for it), and :math:`a`/:math:`r` are whatever absolute/relative gap value is at
hand at the moment of the check. The incumbent-only call site (which
has no real gap value to offer) passes :math:`a = {+\infty}` --
correctly disabling the absolute-gap early-exit -- but currently passes
:math:`r = 1` rather than :math:`+\infty` for the relative side, as a
"gap is certainly not smaller than 100%" stand-in. This is a *latent*
edge case, not a bug in ordinary use: since ``mip:plateauwarmuprelgap``
is a fraction and is realistically always set well below `1`, the
sentinel is never actually reached in practice -- but a
``mip:plateauwarmuprelgap`` set to `1.0` or higher would incorrectly
end warmup on the very first incumbent report, before any gap is
actually known. Using :math:`+\infty` for both parameters at this call
site would close that edge case entirely.

The time-based part exists to avoid judging the search before it has
had a fair chance to get started: without it, a solve that takes
:math:`T_{\text{plateau}}` seconds just to find its *first* incumbent
would be indistinguishable from one that plateaued immediately. The
gap-based part is a convenience on top: if the search already starts
out reasonably close to optimal, there's no need to wait out the full
time-based grace period before plateau detection kicks in. The
warm-up clock :math:`t_0` starts when the feature is armed (i.e. at
solve start, alongside the interrupter setup), not when the first
qualifying report arrives -- so, unlike :math:`t_{\text{last}}`, it is
unconditional on any solver activity and only needs a single fixed
reference point.


.. _plateauPeriodicUtility:

The periodic-callback utility
================================

Not every native callback context has a fresh incumbent or gap value
to report -- e.g. Gurobi's ``GRB_CB_POLLING`` is a generic, always-on
context fired regardless of algorithm phase, with no MIP-specific data
attached at all. For contexts like this, drivers can call
``CheckTimeoutForPlateau()`` instead of ``ReportIncumbentForPlateau()``/
``ReportGapForPlateau()``: it evaluates the same shared stopping
condition (including the warmup gate) without requiring a fresh value. 
This lets a driver keep the plateau clock being checked even during long 
stretches with no new incumbent, rather than only ever re-evaluating the stop
condition as a side effect of ``ReportIncumbentForPlateau()``/
``ReportGapForPlateau()``.
