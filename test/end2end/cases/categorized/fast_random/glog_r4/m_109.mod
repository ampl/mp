var x{1..3} integer >= -2 <= 1;

minimize o: -4*x[1] - 5*x[2] + 5*x[3];

s.t. c1: alldiff(x[2], x[1]);
