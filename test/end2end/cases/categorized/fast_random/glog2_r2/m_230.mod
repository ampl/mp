var x{1..4} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (5)*x[2] + (5)*x[3] + (2)*x[4];
s.t. c1: (x[3] != 3) and (x[2] != 0) ==> (x[1] >= 3);
s.t. c2: atleast 1 {i in 1..4} (x[i] >= 2);
s.t. c3: abs(x[2] - x[3]) != 2;
s.t. c4: (if x[2] >= 1 then x[4] else x[2]) > 3;
