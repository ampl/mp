var x{1..4} integer >= 0 <= 3;

maximize o: 2*x[1] + 6*x[2] + 6*x[3] + 5*x[4];

s.t. c1: abs(x[2] - x[4]) > 3;
s.t. c2: ((if x[4] >= 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[3] == 0) and (x[3] != 3) ==> (x[3] >= 0);
s.t. c4: x[4] == 2;
s.t. c5: x[3] == 0;
