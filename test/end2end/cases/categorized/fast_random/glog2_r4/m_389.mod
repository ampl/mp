var x{1..5} integer >= -2 <= 2;

minimize o: -2*x[1] - 2*x[2] + 0*x[3] + x[4] + 4*x[5];

s.t. c1: alldiff(x[1], x[2], x[4]);
s.t. c2: max(x[5], x[3]) - min(x[2], x[2]) < 1;
s.t. c3: alldiff(x[1], x[2], x[4]);
s.t. c4: max(x[5], x[3]) - min(x[2], x[2]) < 1;
