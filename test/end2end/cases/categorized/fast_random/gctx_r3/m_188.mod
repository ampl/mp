var x{1..3} integer >= 0 <= 2;
minimize o: (2)*x[1] + (1)*x[2] + (-1)*x[3];
s.t. c1: (if max(x[1], x[2]) > 0 then 1 else 0) + (if x[3] <= 1 then 1 else 0) + (if max(x[1], x[2]) > 0 then 1 else 0) >= 1;
s.t. c2: (x[3] <= 1) ==> ((max(x[1], x[2]) > 0) or (max(x[1], x[2]) > 0));
s.t. c3: (max(x[1], x[2]) > 0) ==> (x[3] <= 1);
