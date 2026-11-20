var x{1..4} integer >= 0 <= 4;

minimize o: -2*x[1] - 2*x[2] - 5*x[3] + 6*x[4];

s.t. c1: not ((x[3] <= 2) and (x[1] <= 0));
s.t. c2: ((if x[4] != 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
s.t. c3: ((if x[2] < 3 then 1 else 0) + (if x[4] >= 3 then 1 else 0)) >= 1;
s.t. c4: (x[1] == 0) or (x[1] != 2) or (x[1] != 1);
s.t. c5: (x[3] < 3) or (x[2] >= 1) or (x[2] > 4);
