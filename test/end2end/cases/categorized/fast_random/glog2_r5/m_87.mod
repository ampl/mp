var x{1..4} integer >= 0 <= 2;
maximize o: (34)*x[1] + (-2000)*x[2] + (17)*x[3] + (34)*x[4];
s.t. c1: (x[2] < 0) and (x[2] > 2) ==> (x[3] >= 0);
s.t. c2: ((if x[2] == 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c3: (if x[4] >= 2 then x[2] else x[4]) > 2;
