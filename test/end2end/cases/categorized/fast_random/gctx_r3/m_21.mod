var x{1..3} integer >= 0 <= 2;

minimize o: 3*x[1] - 4*x[2] - 4*x[3];

s.t. c1: (abs(x[2] - x[1]) <= 1) ==> (x[2] < x[3]);
s.t. c2: not (x[2] < x[3]);
s.t. c3: (x[1] == 0) and ((x[2] <= 0) or not (x[2] < x[3]));
