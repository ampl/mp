var x{1..4} integer >= 0 <= 4;

maximize o: -2*x[1] + 0*x[2] + 0*x[3] + 5*x[4];

s.t. c1: abs(x[1] - x[2]) != 1;
s.t. c2: ((if x[4] == 3 then 1 else 0) + (if x[4] >= 3 then 1 else 0)) >= 1;
s.t. c3: ((x[4] > 0) ==> (x[2] != 4)) and ((x[2] != 4) ==> (x[4] >= 0));
