var x{1..5} integer >= -2 <= 2;

maximize o: 0*x[1] - 4*x[2] + 6*x[3] + x[4] - 2*x[5];

s.t. c1: not ((x[4] < 0) and (x[5] < -2));
s.t. c2: (x[5] >= -2) or (x[2] == -1) or (x[4] > 2);
s.t. c3: (x[5] > 0) ==> (x[4] <= -1) else (x[4] < 1);
