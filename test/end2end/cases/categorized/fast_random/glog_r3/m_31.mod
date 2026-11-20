var x{1..5} integer >= 0 <= 4;

minimize o: x[1] - 4*x[2] - 2*x[3] + 4*x[4] + 2*x[5];

s.t. c1: (x[1] >= 4) <==> (x[3] + x[2] == 3);
s.t. c2: max(x[5], x[2]) - min(x[5], x[5]) > 2;
s.t. c3: alldiff(x[5], x[2], x[1], x[3]);
s.t. c4: max(x[4], x[3]) - min(x[2], x[1]) == 1;
