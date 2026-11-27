var x{1..4} integer >= 0 <= 2;

maximize o: -4*x[1] + 0*x[2] + 0*x[3] + 4*x[4] + 6*(if x[3] > 1 then 1 else 0);

s.t. c1: not (x[3] != 2);
s.t. c2: (x[3] != 2) or (min(x[4], x[1]) >= 1);
s.t. c3: not ((x[3] != 2) and (min(x[4], x[1]) >= 1));
s.t. c4: (x[3] != 2) or (x[3] != 2);
