var x{1..5} integer >= 0 <= 1;

minimize o: 2*x[1] - x[2] + 3*x[3] - 2*x[4] - 4*x[5];

s.t. c1: atleast 1 {i in 1..5} (x[i] >= 1);
s.t. c2: not ((x[2] > 0) and (x[5] > 1));
s.t. c3: (x[4] >= 0) <==> (x[4] + x[4] == 1);
s.t. c4: alldiff(x[3], x[2], x[4]);
