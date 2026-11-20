var x{1..3} integer >= 0 <= 3;

maximize o: 0*x[1] + 3*x[2] + 0*x[3];

s.t. c1: (x[3] + x[3] <= 1) ==> ((x[2] != 0) or (x[3] + x[3] <= 1));
s.t. c2: (x[2] != 0) and ((x[3] + x[3] <= 1) or not (x[2] != 0));
