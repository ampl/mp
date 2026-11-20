var x{1..5} integer >= 0 <= 2;

minimize o: -5*x[1] - 4*x[2] + 2*x[3] + 4*x[4] + 3*x[5];

s.t. c1: ((if x[3] > 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c2: atleast 3 {i in 1..5} (x[i] >= 2);
