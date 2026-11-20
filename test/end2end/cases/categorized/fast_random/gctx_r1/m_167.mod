var x{1..4} integer >= 0 <= 3;

minimize o: 4*x[1] + 3*x[2] - 5*x[3] - 2*x[4] - 5*(if x[3] == 2 then 1 else 0);

s.t. c1: (x[3] == 2) and ((x[2] > 3) or not (x[3] == 2));
s.t. c2: not (x[2] > 3);
