var x{1..5} integer >= 0 <= 3;

maximize o: -x[1] + x[2] - 2*x[3] - x[4] + 4*x[5];

s.t. c1: (x[1] >= 0) ==> (x[4] > x[5]);
s.t. c2:
   (if x[4] > x[5] then x[3] else 0) + (if x[4] > x[5] then x[5] else 0) >= 1;
