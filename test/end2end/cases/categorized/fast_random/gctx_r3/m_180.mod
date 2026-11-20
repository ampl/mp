var x{1..3} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] - x[3];

s.t. c1: (x[3] < 2) ==> (x[3] != x[2]) else (x[3] != x[2]);
s.t. c2:
   (if x[3] < 2 then 1 else 0) + (if x[3] < 2 then 1 else 0)
   + (if x[3] < 2 then 1 else 0)
   >= 3;
