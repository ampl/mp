var x{1..3} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (6)*x[2] + (-1)*x[3];
s.t. c1: ((if x[2] >= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c2: atleast 2 {i in 1..3} (x[i] >= 2);
