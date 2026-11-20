var x{1..6} integer >= 0 <= 3;

maximize o: 5*x[1] + 6*x[2] - 2*x[3] - x[4] + 4*x[5] + 4*x[6];

s.t. c1: (x[1] < 3) ==> (x[5] > 3);
s.t. c2: not ((x[4] < 1) and (x[4] <= 2));
s.t. c3: abs(x[2] - x[2]) > 1;
