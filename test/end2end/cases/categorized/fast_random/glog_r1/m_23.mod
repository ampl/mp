var x{1..4} integer >= 0 <= 4;

maximize o: -4*x[1] + 6*x[2] - 3*x[3] + x[4];

s.t. c1: (x[2] > 1) and (x[3] <= 4) ==> (x[4] >= 1);
s.t. c2: not ((x[2] == 2) and (x[1] == 2));
s.t. c3: (x[4] < 1) or (x[2] > 4) or (x[1] <= 3);
