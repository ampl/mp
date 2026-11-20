var x{1..6} integer >= 0 <= 3;

minimize o: -4*x[1] + x[2] - 2*x[3] - 3*x[4] - 3*x[5] - x[6];

s.t. c1: alldiff(x[1], x[4], x[2]);
s.t. c2: (x[3] != 1) or (x[5] != 0) or (x[1] != 2);
s.t. c3: ((x[6] < 3) ==> (x[2] > 1)) and ((x[2] > 1) ==> (x[6] >= 3));
s.t. c4: alldiff(x[2], x[4], x[1], x[3]);
