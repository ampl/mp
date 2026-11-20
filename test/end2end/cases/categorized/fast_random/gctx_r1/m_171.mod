var x{1..5} integer >= -1 <= 1;

maximize o: 0*x[1] - x[2] + 4*x[3] + 2*x[4] + 6*x[5];

s.t. c1: (x[3] + x[1] < 0) or (x[1] < -1);
s.t. c2:
   (if x[1] < -1 then x[3] else 0) + (if x[3] + x[1] < 0 then x[5] else 0) < -1;
s.t. c3: (x[3] + x[1] < 0) ==> (x[1] < -1) else (x[3] + x[1] < 0);
