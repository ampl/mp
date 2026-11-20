var x{1..4} integer >= 0 <= 4;
minimize o: (5)*x[1] + (0)*x[2] + (0)*x[3] + (2)*x[4];
s.t. c1: ((x[4] <= 2) ==> (x[4] != 0)) and ((x[4] != 0) ==> (x[3] >= 2));
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c3: (if x[3] == 1 then x[2] else x[3]) >= 1;
s.t. c4: ((if x[2] == 4 then 1 else 0) + (if x[1] >= 4 then 1 else 0)) >= 1;
s.t. c5: ((x[4] <= 2) ==> (x[4] != 0)) and ((x[4] != 0) ==> (x[3] >= 2));
s.t. c6: ((if x[2] == 4 then 1 else 0) + (if x[1] >= 4 then 1 else 0)) >= 1;
