var x{1..4} integer >= 0 <= 4;

minimize o: 3*x[1] - 2*x[2] + 0*x[3] + 6*x[4];

s.t. c1: max(x[2], x[1]) - min(x[4], x[2]) <= 2;
s.t. c2: (x[3] == 3) and (x[2] >= 3) ==> (x[4] >= 3);
s.t. c3: max(x[3], x[3]) - min(x[3], x[4]) < 0;
