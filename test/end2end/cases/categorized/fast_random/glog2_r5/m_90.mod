var x{1..5} integer >= -2 <= 2;

minimize o: -x[1] + 6*x[2] + 3*x[3] - 4*x[4] - x[5];

s.t. c1: abs(x[4] - x[1]) <= 4;
s.t. c2: not ((x[4] != 2) and (x[1] >= 0));
s.t. c3: (x[2] < -1) ==> (x[5] >= 0) else (x[2] == -2);
s.t. c4: (x[4] == 1) ==> (x[5] <= 1) else (x[1] == 0);
