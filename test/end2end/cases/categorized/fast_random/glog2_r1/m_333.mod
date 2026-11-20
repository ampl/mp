var x{1..6} integer >= 0 <= 4;

maximize o: 4*x[1] + 0*x[2] + 2*x[3] - 5*x[4] + 0*x[5] + x[6];

s.t. c1: not ((x[5] > 1) and (x[1] > 1));
s.t. c2: (x[5] > 1) ==> (x[2] > 1) else (x[6] != 0);
s.t. c3: x[6] == 4;
s.t. c4: x[3] == 0;
