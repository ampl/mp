var x{1..7} integer >= 0 <= 1;

minimize o: -3*x[1] + 0*x[2] + x[3] - 3*x[4] + x[5] - 4*x[6] - x[7];

s.t. c1: (x[6] == 0) <==> (x[7] + x[2] >= 0);
s.t. c2: exists {i in 1..7} x[i] == 0;
s.t. c3: ((x[5] > 0) ==> (x[1] >= 1)) and ((x[1] >= 1) ==> (x[5] >= 0));
s.t. c4: alldiff(x[5], x[7], x[4]);
