var x{1..4} integer >= -2 <= 2;

maximize o: x[1] - 4*x[2] + 6*x[3] + 5*x[4];

s.t. c1: (x[2] == -1) <==> (x[2] + x[2] < -2);
s.t. c2: not ((x[2] != -2) and (x[1] != 2));
