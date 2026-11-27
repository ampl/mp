var x{1..3} integer >= -2 <= 1;

maximize o: -1000*x[1] + 17*x[2] + 51*x[3];

s.t. c1: ((if x[3] == 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
