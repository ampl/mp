var x{1..5} integer >= 0 <= 4;

minimize o: -2*x[1] + 6*x[2] + 4*x[3] - x[4] + 4*x[5];

s.t. c1: alldiff(x[5], x[4]);
s.t. c2: ((if x[4] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c3: alldiff(x[5], x[2]);
s.t. c4: (if x[1] != 1 then x[3] else x[4]) != 3;
s.t. c5: x[2] == 0;
