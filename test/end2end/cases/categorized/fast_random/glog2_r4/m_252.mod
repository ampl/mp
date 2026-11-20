var x{1..4} integer >= 0 <= 2;

maximize o: 4*x[1] + 4*x[2] + 5*x[3] + x[4];

s.t. c1: max(x[2], x[2]) - min(x[2], x[1]) <= 1;
s.t. c2: ((x[1] != 1) ==> (x[1] == 0)) and ((x[1] == 0) ==> (x[1] >= 1));
s.t. c3: ((if x[3] < 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c4: (x[4] == 2) or (x[4] > 1) or (x[4] != 0);
