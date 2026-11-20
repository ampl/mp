var x{1..3} integer >= 0 <= 2;
minimize o: (4)*x[1] + (-3)*x[2] + (4)*x[3];
s.t. c1: (x[1] < 1) and (x[3] != 0) ==> (x[3] >= 1);
s.t. c2: ((if x[2] > 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
