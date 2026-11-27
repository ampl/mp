var x{1..3} integer >= -2 <= 2;

minimize o: 3*x[1] - 3*x[2] + x[3];

s.t. c1: alldiff(x[2], x[3]);
s.t. c2: alldiff(x[1], x[2]);
s.t. c3: alldiff(x[2], x[3]);
s.t. c4: alldiff(x[1], x[2]);
