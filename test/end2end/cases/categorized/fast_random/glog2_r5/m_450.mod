var x{1..4} integer >= -2 <= 2;
maximize o: (2)*x[1] + (-4)*x[2] + (1)*x[3] + (4)*x[4];
s.t. c1: (x[2] < 1) and (x[3] > 0) ==> (x[2] >= 1);
s.t. c2: (x[2] == 2) and (x[1] > -1) ==> (x[2] >= 2);
s.t. c3: (if x[1] <= 2 then x[4] else x[1]) >= -1;
