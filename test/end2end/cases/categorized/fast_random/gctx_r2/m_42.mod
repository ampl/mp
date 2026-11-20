var x{1..3} integer >= 0 <= 2;

minimize o: -4*x[1] + 2*x[2] - 5*x[3];

s.t. c1: (if x[2] >= 0 then x[2] else x[1]) < 2;
s.t. c2: (x[2] >= 0) or (x[2] >= 0);
s.t. c3: (x[2] >= 0) <==> (x[3] <= 0);
