var x{1..3} integer >= 0 <= 2;

maximize o: 5*x[1] - 2*x[2] + 3*x[3];

s.t. c1: (x[1] + x[2] >= 1) <==> (x[1] + x[2] >= 1);
s.t. c2: not ((x[1] != x[1]) and (x[1] != x[1]));
