var x{1..3} integer >= -2 <= 1;
maximize o: (-4)*x[1] + (3)*x[2] + (5)*x[3];
s.t. c1: abs(x[2] - x[2]) <= 1;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (if x[2] == 1 then x[1] else x[3]) != -2;
s.t. c4: ((x[3] == -1) ==> (x[3] <= 0)) and ((x[3] <= 0) ==> (x[1] >= -1));
