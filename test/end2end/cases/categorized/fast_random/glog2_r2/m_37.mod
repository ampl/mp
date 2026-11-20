var x{1..3} integer >= 0 <= 2;

maximize o: 17*x[1] - 999*x[2] - 1000*x[3];

s.t. c1: ((if x[2] >= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
