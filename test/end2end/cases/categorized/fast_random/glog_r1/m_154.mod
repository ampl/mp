var x{1..5} integer >= 0 <= 4;

minimize o: -4*x[1] - 4*x[2] - x[3] - 5*x[4] + x[5];

s.t. c1: (x[4] != 2) and (x[3] > 3) ==> (x[1] >= 2);
s.t. c2: max(x[2], x[2]) - min(x[4], x[2]) <= 4;
s.t. c3: not ((x[5] > 1) and (x[4] <= 2));
