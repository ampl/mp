var x{1..3} integer >= 0 <= 2;

maximize o: -4*x[1] + 4*x[2] - 4*x[3] + 6*(if x[2] != 0 then 1 else 0);

s.t. c1: (x[1] == x[2]) <==> (x[1] == x[2]);
s.t. c2: (x[2] != 0) and ((x[3] <= x[3]) or not (x[3] <= x[3]));
s.t. c3: not (x[2] > 2);
