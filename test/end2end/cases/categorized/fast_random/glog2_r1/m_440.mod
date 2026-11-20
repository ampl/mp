var x{1..5} integer >= -2 <= 2;

maximize o: -4*x[1] + 3*x[2] + 4*x[3] - 4*x[4] + 5*x[5];

s.t. c1: abs(x[4] - x[1]) > 1;
s.t. c2: (x[4] >= -2) and (x[1] != 0) ==> (x[3] >= -2);
