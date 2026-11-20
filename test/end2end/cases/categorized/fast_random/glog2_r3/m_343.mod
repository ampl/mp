var x{1..4} integer >= -2 <= 2;

maximize o: 4*x[1] + 2*x[2] - 3*x[3] + 6*x[4];

s.t. c1: (x[4] == 0) and (x[1] != 1) ==> (x[2] >= 0);
s.t. c2: x[1] == -1;
s.t. c3: x[2] == 0;
