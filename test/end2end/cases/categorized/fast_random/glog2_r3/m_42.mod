var x{1..4} integer >= 0 <= 4;

maximize o: x[1] + 4*x[2] + 5*x[3] - x[4];

s.t. c1: (x[2] < 4) and (x[2] == 1) ==> (x[4] >= 4);
s.t. c2: abs(x[2] - x[2]) > 0;
s.t. c3: ((if x[4] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c4: max(x[2], x[2]) - min(x[2], x[1]) == 3;
