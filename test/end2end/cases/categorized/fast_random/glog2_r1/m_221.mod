var x{1..6} integer >= 0 <= 4;
minimize o: (3)*x[1] + (-2)*x[2] + (-5)*x[3] + (2)*x[4] + (-2)*x[5] + (-1)*x[6];
s.t. c1: ((x[6] >= 0) ==> (x[3] == -1)) and ((x[3] == -1) ==> (x[5] >= 0));
s.t. c2: numberof 5 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 3;
s.t. c3: ((if x[1] != 5 then 1 else 0) + (if x[1] >= 5 then 1 else 0)) >= 1;
