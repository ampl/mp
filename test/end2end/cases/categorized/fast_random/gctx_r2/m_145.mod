var x{1..3} integer >= -1 <= 1;

maximize o: 0*x[1] - x[2] + 6*x[3];

s.t. c1: (x[2] > x[2]) or (x[3] != x[2]);
s.t. c2:
   (if x[2] > x[2] then 1 else 0) + (if x[3] != x[2] then 1 else 0)
   + (if x[2] > x[2] then 1 else 0)
   >= 2;
s.t. c3: (x[2] <= -1) ==> (x[2] <= -1);
