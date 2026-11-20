var x{1..4} integer >= 0 <= 3;

maximize o: 6*x[1] + x[2] + 6*x[3] + 0*x[4];

s.t. c1: (if x[1] != 0 then x[2] else 0) + (if x[1] != 0 then x[2] else 0) <= 0;
s.t. c2: (x[4] + x[1] <= 0) and ((x[4] + x[1] <= 0) or not (x[1] != 0));
s.t. c3:
   (if x[4] + x[1] <= 0 then 1 else 0) + (if x[1] != 0 then 1 else 0)
   + (if x[4] + x[1] <= 0 then 1 else 0)
   >= 1;
