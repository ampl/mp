var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] + 6*x[2] - 3*x[3];

s.t. c1: (abs(x[1] - x[2]) <= 3) ==> (abs(x[2] - x[1]) > 0);
s.t. c2: (x[3] != x[2]) ==> (abs(x[1] - x[1]) < 2);
s.t. c3: not (abs(x[1] - x[1]) < 2);
