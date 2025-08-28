## Test cvt:pre:unnest, cvt:pre:sort

set I default 1..12;

var x {I} >=-20 <= 67 integer;

s.t. ConOr1: (
  (forall {i in I: i mod 4 == 1}
     x[i] >= i+3 && x[i+1] >= i+4 && x[i+2] >= i+5)
     &&
       forall {i in I: i mod 4 == 0} x[i] >= i+3)
       || (
         (forall {i in I: i mod 4 == 0}
           (x[i] >= i+3 && x[i-1] >= i+2))     ## repeat
         &&
           (forall {i in I: i mod 4 == 1}
             x[i] >= i+3 && x[i+2] >= i+5 && x[i+1] >= i+4));

s.t. ConOr21: x[3] <= 7 or x[4] <= 20 or x[1] <= 6;

s.t. ConOr22: (x[3] <= 7 or x[1] <= 6)
               or x[4] <= 20 or (x[1] <= 6 or x[3] <= 7);

minimize O: sum {i in I} x[i];
