var x{1..4} integer >= 0 <= 4;
minimize o: (51)*x[1] + (-1000)*x[2] + (2000)*x[3] + (-333)*x[4];
s.t. c1: ((if x[2] < 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
