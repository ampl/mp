var x{1..6} integer >= 0 <= 4;

minimize o: -4*x[1] - x[2] - x[3] + 4*x[4] + 5*x[5] + 3*x[6];

s.t. c1: alldiff(x[5], x[2]);
s.t. c2: ((x[6] == 1) ==> (x[5] > 1)) and ((x[5] > 1) ==> (x[5] >= 1));
s.t. c3: alldiff(x[6], x[5], x[4]);
s.t. c4: (x[4] == 3) ==> (x[2] <= 1) else (x[6] <= 4);
