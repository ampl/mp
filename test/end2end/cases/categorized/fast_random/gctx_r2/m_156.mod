var x{1..3} integer >= -1 <= 2;

maximize o: 4*x[1] + 0*x[2] - 4*x[3] + 5*(if x[1] == x[1] then 1 else 0);

s.t. c1: (x[1] == x[1]) ==> ((x[1] == x[1]) or (x[1] == x[1]));
s.t. c2: (x[2] <= 2) or (x[1] == x[1]);
s.t. c3: (x[2] <= 2) ==> (x[2] <= 2) else (x[2] <= x[2]);
