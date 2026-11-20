var x{1..3} integer >= 0 <= 3;

minimize o: 3*x[1] + 5*x[2] - 3*x[3] - 5*(if x[2] == 3 then 1 else 0);

s.t. c1: (x[1] != 1) ==> (x[1] != 1);
s.t. c2: not (x[2] == 3);
