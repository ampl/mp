var x{1..3} integer >= 0 <= 2;

maximize o: -4*x[1] - x[2] - 3*x[3];

s.t. c1: (min(x[2], x[1]) > 1) or (x[2] + x[2] > 1);
s.t. c2: (x[2] + x[2] > 1) ==> (min(x[2], x[1]) > 1) else (min(x[2], x[1]) > 1);
s.t. c3: (x[2] + x[2] > 1) or (min(x[2], x[1]) > 1);
