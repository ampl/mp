var x{1..4} integer >= 0 <= 2;

minimize o: -4*x[1] - x[2] + 6*x[3] - 2*x[4];

s.t. c1: (x[4] >= 2) or (x[2] == 1) or (x[2] != 2);
s.t. c2: max(x[3], x[1]) - min(x[3], x[3]) == 0;
s.t. c3: (x[2] == 0) and (x[4] < 0) ==> (x[1] >= 0);
