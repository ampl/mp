var x{1..3} integer >= 0 <= 4;

minimize o: 3*x[1] - 5*x[2] + 6*x[3];

s.t. c1: abs(x[3] - x[1]) == 2;
s.t. c2: not ((x[3] == 3) and (x[1] == 0));
s.t. c3: alldiff(x[1], x[3], x[2]);
s.t. c4: ((if x[3] >= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c5: (x[1] <= 1) or (x[2] < 3) or (x[3] == 4);
s.t. c6: not ((x[3] < 1) and (x[3] != 0));
