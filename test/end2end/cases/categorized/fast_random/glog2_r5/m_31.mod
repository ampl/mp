var x{1..5} integer >= 0 <= 4;

maximize o: -2*x[1] - 3*x[2] + x[3] - x[4] - 4*x[5];

s.t. c1: ((if x[4] > 5 then 1 else 0) + (if x[2] >= 5 then 1 else 0)) >= 1;
s.t. c2: ((x[1] > 5) ==> (x[5] > 5)) and ((x[5] > 5) ==> (x[4] >= 5));
