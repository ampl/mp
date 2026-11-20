var x{1..5} integer >= 0 <= 4;

minimize o: 4*x[1] - 2*x[2] - x[3] + 0*x[4] - 2*x[5];

s.t. c1: (x[2] != 1) and (x[4] > 0) ==> (x[3] >= 1);
s.t. c2: max(x[3], x[4]) - min(x[3], x[3]) == 2;
s.t. c3: not ((x[2] < 0) and (x[4] != 3));
s.t. c4: x[3] == 1;
