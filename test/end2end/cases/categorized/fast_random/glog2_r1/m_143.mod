var x{1..4} integer >= 0 <= 2;

minimize o: 6*x[1] - 3*x[2] + x[3] + 6*x[4];

s.t. c1: (x[3] > 1) and (x[4] > 2) ==> (x[2] >= 1);
s.t. c2: not ((x[3] != 2) and (x[4] == 1));
s.t. c3: max(x[2], x[4]) - min(x[2], x[1]) < 1;
s.t. c4: (x[3] != 0) ==> (x[4] <= 1);
s.t. c5: x[4] == 2;
s.t. c6: x[2] == 1;
