var x{1..6} integer >= 0 <= 1;

minimize o: 6*x[1] + 0*x[2] - 3*x[3] + 2*x[4] - 2*x[5] + 0*x[6];

s.t. c1: ((x[3] != 1) ==> (x[4] > 0)) and ((x[4] > 0) ==> (x[3] >= 1));
s.t. c2: abs(x[1] - x[2]) < 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 1;
