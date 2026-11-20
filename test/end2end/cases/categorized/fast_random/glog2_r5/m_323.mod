var x{1..6} integer >= 0 <= 3;

maximize o: 0*x[1] + 3*x[2] + 5*x[3] + 4*x[4] - 2*x[5] + 2*x[6];

s.t. c1: (x[6] <= 2) or (x[5] == 2) or (x[4] > 1);
s.t. c2: not ((x[6] == 3) and (x[4] >= 0));
s.t. c3: x[5] == 3;
