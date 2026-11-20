var x{1..6} integer >= 0 <= 3;

maximize o: -5*x[1] + 0*x[2] + 4*x[3] - 2*x[4] + 0*x[5] - 2*x[6];

s.t. c1: (x[3] >= 1) ==> (x[2] > 2);
s.t. c2: not ((x[5] != 3) and (x[1] == 3));
s.t. c3: ((if x[2] <= 3 then 1 else 0) + (if x[6] >= 3 then 1 else 0)) >= 1;
s.t. c4: (x[2] < 2) ==> (x[3] < 1) else (x[1] == 0);
s.t. c5: (x[3] >= 1) ==> (x[2] > 2);
s.t. c6: (x[2] < 2) ==> (x[3] < 1) else (x[1] == 0);
