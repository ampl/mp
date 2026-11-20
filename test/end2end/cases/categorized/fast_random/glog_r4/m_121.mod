var x{1..6} integer >= 0 <= 4;

maximize o: -5*x[1] - 4*x[2] + 5*x[3] + 5*x[4] + 2*x[5] - 4*x[6];

s.t. c1: ((x[4] == 2) ==> (x[3] != 3)) and ((x[3] != 3) ==> (x[2] >= 2));
s.t. c2: ((if x[4] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c3: ((if x[2] <= 3 then 1 else 0) + (if x[4] >= 3 then 1 else 0)) >= 1;
