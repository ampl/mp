var x{1..4} integer >= 0 <= 4;

minimize o: x[1] - 5*x[2] + 5*x[3] - x[4];

s.t. c1: (x[3] <= 1) <==> (x[3] + x[1] == 2);
s.t. c2: alldiff(x[3], x[4]);
