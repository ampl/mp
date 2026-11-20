var x{1..4} integer >= -1 <= 2;
maximize o: (2)*x[1] + (-2)*x[2] + (5)*x[3] + (-5)*x[4];
s.t. c1: (if x[3] == 1 then x[2] else 0) + (if x[3] == 1 then x[1] else 0) < -1;
s.t. c2: (x[3] == 1) ==> ((x[3] == 1) or (x[3] == 1));
s.t. c3: (x[3] == 1) ==> (x[1] <= x[4]);
s.t. c4: (x[3] == 1) and ((x[1] <= x[4]) or not (x[3] == 1));
