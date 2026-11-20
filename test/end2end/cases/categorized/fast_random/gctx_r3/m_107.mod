var x{1..3} integer >= 0 <= 2;

maximize o: -x[1] + 0*x[2] + 4*x[3] - 5*(if x[1] != 0 then 1 else 0);

s.t. c1: (x[1] + x[2] == 2) or (x[1] + x[2] == 2);
s.t. c2: (x[1] != 0) ==> ((x[1] != 0) or (x[1] > 0));
