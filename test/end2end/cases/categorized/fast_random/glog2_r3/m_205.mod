var x{1..3} integer >= 0 <= 3;
maximize o: (5)*x[1] + (6)*x[2] + (4)*x[3];
s.t. c1: (x[2] > 2) or (x[3] >= 3) or (x[2] <= 2);
s.t. c2: not ((x[1] >= 0) and (x[1] < 3));
s.t. c3: (if x[3] >= 2 then x[1] else x[1]) != 2;
s.t. c4: not ((x[1] >= 3) and (x[2] != 3));
s.t. c5: ((x[2] != 3) ==> (x[1] != 3)) and ((x[1] != 3) ==> (x[1] >= 3));
