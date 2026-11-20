var x{1..4} integer >= -1 <= 1;

maximize o: 5*x[1] + 4*x[2] - 4*x[3] + 4*x[4];

s.t. c1: (abs(x[3] - x[1]) > 2) <==> (min(x[4], x[3]) != -1);
s.t. c2: (abs(x[3] - x[1]) > 2) <==> (min(x[4], x[3]) != -1);
