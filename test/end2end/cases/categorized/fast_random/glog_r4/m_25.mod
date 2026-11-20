var x{1..4} integer >= -2 <= 2;
minimize o: (-4)*x[1] + (-1)*x[2] + (1)*x[3] + (-2)*x[4];
s.t. c1: ((if x[4] != 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[3] >= 2) ==> (x[4] == 2);
