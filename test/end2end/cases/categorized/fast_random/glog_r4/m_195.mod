var x{1..6} integer >= 0 <= 3;

minimize o: 2*x[1] + 6*x[2] + 2*x[3] - 4*x[4] - 3*x[5] + 2*x[6];

s.t. c1: atleast 1 {i in 1..6} (x[i] >= 0);
s.t. c2: not ((x[2] > 3) and (x[4] <= 3));
s.t. c3: ((x[1] <= 2) ==> (x[3] == 0)) and ((x[3] == 0) ==> (x[6] >= 2));
s.t. c4: alldiff(x[4], x[5], x[2]);
