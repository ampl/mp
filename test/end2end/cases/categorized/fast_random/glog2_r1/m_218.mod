var x{1..5} integer >= 0 <= 1;

minimize o: 2*x[1] + 4*x[2] - 4*x[3] + 5*x[4] + 6*x[5];

s.t. c1: (x[2] >= 0) ==> (x[3] < 0) else (x[1] == 1);
s.t. c2: max(x[4], x[5]) - min(x[5], x[5]) > 0;
s.t. c3: alldiff(x[1], x[3], x[4]);
s.t. c4: ((x[4] < 0) ==> (x[5] > 1)) and ((x[5] > 1) ==> (x[2] >= 0));
