var x{1..4} integer >= -2 <= 0;

minimize o: -2*x[1] + x[2] + x[3] + 5*x[4];

s.t. c1: ((if x[4] < 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: atleast 2 {i in 1..4} (x[i] >= -2);
