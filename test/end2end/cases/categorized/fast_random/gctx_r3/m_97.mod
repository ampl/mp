var x{1..5} integer >= -1 <= 1;

maximize o: -2*x[1] - x[2] + 0*x[3] + 2*x[4] + 2*x[5];

s.t. c1: (x[4] < 0) ==> (x[3] < 0);
s.t. c2:
   (if x[4] < 0 then x[5] else 0) + (if x[1] >= x[4] then x[5] else 0) <= 1;
