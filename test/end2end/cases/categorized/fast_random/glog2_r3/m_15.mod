var x{1..4} integer >= 0 <= 4;
maximize o: (3)*x[1] + (3)*x[2] + (0)*x[3] + (5)*x[4];
s.t. c1: (x[2] != 3) or (x[2] != 2) or (x[1] < 3);
s.t. c2: (x[2] <= 0) or (x[3] <= 1) or (x[1] != 2);
s.t. c3: (if x[2] != 1 then x[3] else x[1]) != 1;
s.t. c4: not ((x[1] > 1) and (x[1] < 4));
s.t. c5: ((if x[3] > 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
s.t. c6: (x[3] <= 4) and (x[2] >= 2) ==> (x[1] >= 4);
s.t. c7: ((x[3] < 0) ==> (x[3] == 2)) and ((x[3] == 2) ==> (x[4] >= 0));
