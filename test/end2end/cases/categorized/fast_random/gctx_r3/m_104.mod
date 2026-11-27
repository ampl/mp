var x{1..4} integer >= -1 <= 2;

maximize o: 4*x[1] - 5*x[2] + 3*x[3] + 0*x[4];

s.t. c1: (x[1] < -1) <==> (x[4] >= 2);
s.t. c2:
   (if x[1] < -1 then x[1] else 0) + (if x[3] + x[1] <= -1 then x[1] else 0)
   >= 0;
