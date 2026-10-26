# -------------------------------------------------------------
# int_ops_02.mod
# floor, ceil, round (ties to even), trunc, round(x, n), and
# div/mod of real numbers, on continuous variables.
# Each term attains its maximum at a closed end of a piece:
#   floor(a) - 0.01a:      a = 4         ->  3.96
#   ceil(b) + 0.01b:       b = 3.3       ->  4.033
#   round(c) + 0.01c:      c = 2.6       ->  3.026
#   trunc(d) - 0.01d:      d = 2         ->  1.98
#   e div 2.5 - e mod 2.5 - 0.01e:  e = 7.5 -> 2.925
#   round(g, 1) + 0.01g:   g = 0.33      ->  0.3033
#     (decimal ties are not exact in binary: avoided)
#   round(h) - 0.01h:      h = 1.5 (ties to even: 2) -> 1.985
# -------------------------------------------------------------

var a >= -3.7, <= 4.2;
var b >= -2.5, <= 3.3;
var c >= -1.2, <= 2.6;
var d >= -3.3, <= 2.2;
var e >= 1,    <= 9.5;
var g >= 0.12, <= 0.33;
var h >= -0.5, <= 2.5;

maximize Obj:
    floor(a) - 0.01*a
  + ceil(b) + 0.01*b
  + round(c) + 0.01*c
  + trunc(d) - 0.01*d
  + (e div 2.5) - (e mod 2.5) - 0.01*e
  + round(g, 1) + 0.01*g
  + round(h) - 0.01*h;
