# -------------------------------------------------------------
# chk_logtol_01.mod
# sol:chk:logtol: x is fixed by a constraint 5e-8 below the
# integer 3. A solver may return floor(x) = 3, within its
# feasibility tolerance. The solution check evaluates floor() of
# values within sol:chk:logtol (default 1e-6) of an integer at
# that integer: no violation. With sol:chk:logtol=0 (exact, as
# AMPL): floor(x) = 2, a violation.
# (The check of a comparison itself measures its violation with
# sol:chk:feastol; logtol applies to truth values computed for
# expressions, such as floor, ceil, round and numberof.)
# -------------------------------------------------------------
var x >= 0 <= 10;

maximize o: floor(x);

s.t. fx: x = 3 - 5e-8;
