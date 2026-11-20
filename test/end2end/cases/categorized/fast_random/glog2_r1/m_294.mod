var x{1..5} integer >= -2 <= 3;

maximize o: x[1] + 0*x[2] - 4*x[3] + 3*x[4] + 6*x[5];

s.t. c1: (x[4] > 2) and (x[5] != 0) ==> (x[5] >= 2);
s.t. c2: abs(x[3] - x[2]) > 5;
