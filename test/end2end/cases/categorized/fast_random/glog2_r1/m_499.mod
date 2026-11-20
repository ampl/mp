var x{1..6} integer >= 0 <= 4;

maximize o: 3*x[1] + 5*x[2] - x[3] + 3*x[4] - 4*x[5] + 5*x[6];

s.t. c1: (x[1] <= 3) ==> (x[1] > 0) else (x[1] < 3);
s.t. c2: ((if x[2] == 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
s.t. c3: (x[1] <= 3) ==> (x[1] > 0) else (x[1] < 3);
s.t. c4: ((if x[2] == 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
