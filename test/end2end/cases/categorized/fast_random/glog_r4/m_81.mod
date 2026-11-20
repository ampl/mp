var x{1..4} integer >= 0 <= 2;
minimize o: (-2)*x[1] + (-5)*x[2] + (-3)*x[3] + (1)*x[4];
s.t. c1: ((x[1] > 1) ==> (x[4] >= 2)) and ((x[4] >= 2) ==> (x[2] >= 1));
s.t. c2: ((if x[1] != 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: abs(x[3] - x[1]) > 1;
s.t. c4: (x[4] >= 0) and (x[1] < 2) ==> (x[1] >= 0);
