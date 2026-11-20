var x{1..6} integer >= 0 <= 3;

minimize o: 5*x[1] - x[2] - 4*x[3] - 4*x[4] - 2*x[5] + 0*x[6];

s.t. c1: abs(x[2] - x[5]) == 0;
s.t. c2: alldiff(x[6], x[5]);
s.t. c3: (x[5] > 0) and (x[6] <= 3) ==> (x[3] >= 0);
