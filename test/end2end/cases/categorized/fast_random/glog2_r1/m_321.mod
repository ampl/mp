var x{1..4} integer >= 0 <= 2;

maximize o: 5*x[1] - 3*x[2] + 4*x[3] + 4*x[4];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4]) < 3;
s.t. c2: abs(x[1] - x[3]) >= 1;
s.t. c3: (x[1] < -1) ==> (x[2] != 3);
s.t. c4: ((if x[3] != 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
