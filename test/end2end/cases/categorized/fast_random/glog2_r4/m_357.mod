var x{1..5} integer >= -2 <= 1;

maximize o: -3000*x[1] - 1000*x[2] + 17*x[3] - 1000*x[4] + 51*x[5];

s.t. c1: ((if x[2] < -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
