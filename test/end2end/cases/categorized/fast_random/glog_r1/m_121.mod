var x{1..3} integer >= 0 <= 3;

maximize o: -x[1] - 2*x[2] + 2*x[3];

s.t. c1: (x[1] <= 2) and (x[2] >= 0) ==> (x[1] >= 2);
