var x{1..5} integer >= 0 <= 4;

maximize o: -2*x[1] + 0*x[2] - 4*x[3] - 3*x[4] - 5*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) == 2;
s.t. c2: abs(x[3] - x[3]) > 4;
