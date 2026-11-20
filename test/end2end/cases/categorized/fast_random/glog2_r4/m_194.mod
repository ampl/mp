var x{1..5} integer >= -2 <= 3;
maximize o: (5)*x[1] + (-2)*x[2] + (5)*x[3] + (0)*x[4] + (3)*x[5];
s.t. c1: not ((x[1] > -1) and (x[3] >= 0));
s.t. c2: ((x[1] != 2) ==> (x[5] > 2)) and ((x[5] > 2) ==> (x[5] >= 2));
s.t. c3: ((if x[4] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c4: (if x[1] <= 3 then x[5] else x[5]) > 1;
