var x{1..5} integer >= 0 <= 3;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: ((if x[5] != 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
s.t. c2: ((x[3] != 1) ==> (x[5] < 2)) and ((x[5] < 2) ==> (x[5] >= 1));
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) == 3;
