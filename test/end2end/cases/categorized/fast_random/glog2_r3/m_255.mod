var x{1..4} integer >= -2 <= 1;

maximize o: -4*x[1] - 3*x[2] - 4*x[3] + x[4];

s.t. c1: (x[2] == -2) or (x[3] < -2) or (x[2] <= 0);
s.t. c2: (x[4] > -1) ==> (x[2] >= 1) else (x[4] != 0);
s.t. c3: (x[1] < 0) or (x[2] >= 1) or (x[3] <= 0);
s.t. c4: ((if x[4] == -1 then 1 else 0) + (if x[4] >= -1 then 1 else 0)) >= 1;
s.t. c5: ((if x[1] > 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c6: not ((x[2] <= 0) and (x[4] < -1));
