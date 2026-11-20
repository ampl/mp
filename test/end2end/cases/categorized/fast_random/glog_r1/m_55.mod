var x{1..5} integer >= 0 <= 2;

minimize o: x[1] - 2*x[2] - x[3] - 3*x[4] + 4*x[5];

s.t. c1: (x[4] != 1) ==> (x[4] != 0) else (x[3] < 0);
s.t. c2: (x[1] > 1) and (x[3] < 1) ==> (x[2] >= 1);
s.t. c3: abs(x[4] - x[3]) < 1;
