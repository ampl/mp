var x{1..4} integer >= -2 <= 2;

maximize o: 3*x[1] - x[2] - 2*x[3] + 4*x[4];

s.t. c1: alldiff(x[4], x[2]);
s.t. c2: (x[1] < -2) <==> (x[3] + x[4] != 1);
