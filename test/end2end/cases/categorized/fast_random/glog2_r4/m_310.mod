var x{1..6} integer >= -2 <= 0;

minimize o: 5*x[1] - 3*x[2] + x[3] + 0*x[4] + 5*x[5] + 4*x[6];

s.t. c1: alldiff(x[1], x[6]);
s.t. c2: (x[4] < 0) or (x[6] < -1) or (x[3] != 0);
s.t. c3: (x[1] < -1) <==> (x[6] + x[5] >= 0);
