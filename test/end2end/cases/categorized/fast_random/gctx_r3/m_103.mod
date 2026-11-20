var x{1..4} integer >= 0 <= 3;

maximize o: x[1] - x[2] + 0*x[3] + 0*x[4];

s.t. c1: (min(x[1], x[2]) > 1) ==> (x[2] >= x[3]);
s.t. c2: (if min(x[1], x[2]) > 1 then x[1] else x[2]) > 2;
s.t. c3: (min(x[1], x[2]) > 1) or (min(x[1], x[2]) > 1);
s.t. c4:
   (if x[1] > x[2] then x[1] else 0) + (if x[2] >= x[3] then x[3] else 0) < 2;
