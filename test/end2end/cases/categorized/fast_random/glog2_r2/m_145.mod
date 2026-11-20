var x{1..3} integer >= 0 <= 2;
maximize o: (3)*x[1] + (5)*x[2] + (2)*x[3];
s.t. c1: not ((x[1] != 1) and (x[2] > 0));
s.t. c2: ((x[3] != 2) ==> (x[2] > 0)) and ((x[2] > 0) ==> (x[2] >= 2));
s.t. c3: (if x[1] <= 0 then x[1] else x[1]) == 1;
s.t. c4: ((if x[2] > 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c5: abs(x[1] - x[3]) == 2;
