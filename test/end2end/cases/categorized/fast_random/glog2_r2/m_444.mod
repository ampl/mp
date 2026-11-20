var x{1..5} integer >= -3 <= 2;

maximize o: 3*x[1] + 5*x[2] + 6*x[3] + x[4] - 2*x[5];

s.t. c1: abs(x[5] - x[1]) >= 4;
s.t. c2: (x[5] != 1) ==> (x[1] >= -3);
