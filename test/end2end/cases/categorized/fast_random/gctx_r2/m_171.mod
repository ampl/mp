var x{1..3} integer >= -1 <= 2;

maximize o: 5*x[1] - 3*x[2] + 4*x[3] - 2*(if x[3] + x[3] >= 2 then 1 else 0);

s.t. c1: not ((x[2] >= 1) and (x[3] + x[3] >= 2));
s.t. c2: (x[2] <= x[2]) and ((x[2] <= x[2]) or not (x[2] <= x[2]));
s.t. c3: not (x[2] >= 1);
s.t. c4: (x[2] >= 1) ==> (x[3] + x[3] >= 2) else (x[2] <= x[2]);
