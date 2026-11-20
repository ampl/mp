var x{1..3} integer >= 0 <= 2;

maximize o: -2*x[1] - 2*x[2] - 5*x[3];

s.t. c1: ((if x[1] == 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
