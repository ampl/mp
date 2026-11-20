var x{1..7} integer >= 0 <= 1;

maximize o: 3*x[1] + 3*x[2] + 0*x[3] - 3*x[4] + 2*x[5] + 0*x[6] - x[7];

s.t. c1: ((if x[5] <= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c2: not ((x[1] >= 0) and (x[5] < 1));
s.t. c3: (x[3] >= 1) and (x[5] != 1) ==> (x[6] >= 1);
s.t. c4: ((x[7] >= 1) ==> (x[2] != 0)) and ((x[2] != 0) ==> (x[6] >= 1));
