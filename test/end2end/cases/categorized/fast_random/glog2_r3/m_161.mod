var x{1..6} integer >= -2 <= 0;

minimize o: 3*x[1] + 5*x[2] + 2*x[3] + 2*x[4] - 5*x[5] - 4*x[6];

s.t. c1: atleast 2 {i in 1..6} (x[i] >= -3);
s.t. c2: alldiff(x[6], x[2]);
