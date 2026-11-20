var x{1..3} integer >= 0 <= 4;

maximize o: 2*x[1] + 2*x[2] + 0*x[3];

s.t. c1: abs(x[3] - x[2]) != 2;
s.t. c2: alldiff(x[1], x[3]);
s.t. c3: abs(x[3] - x[2]) != 2;
s.t. c4: alldiff(x[1], x[3]);
