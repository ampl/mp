var x{1..3} integer >= 0 <= 3;

minimize o: -2*x[1] - 3*x[2] - 4*x[3];

s.t. c1: max(x[3], x[3]) - min(x[3], x[3]) != 0;
s.t. c2: abs(x[2] - x[2]) != 1;
s.t. c3: (x[1] == 3) ==> (x[1] >= 1) else (x[2] > 1);
