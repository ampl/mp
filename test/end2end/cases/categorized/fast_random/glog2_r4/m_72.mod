var x{1..5} integer >= 0 <= 4;

maximize o: -4*x[1] + 3*x[2] - x[3] + 5*x[4] + 0*x[5];

s.t. c1: alldiff(x[3], x[1], x[4], x[2]);
s.t. c2: atleast 1 {i in 1..5} (x[i] >= 4);
