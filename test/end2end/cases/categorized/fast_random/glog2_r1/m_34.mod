var x{1..5} integer >= -2 <= 2;

maximize o: -x[1] - 3*x[2] - 5*x[3] + x[4] - 2*x[5];

s.t. c1: ((if x[4] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
