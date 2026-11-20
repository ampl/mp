var x{1..4} integer >= 0 <= 3;
maximize o: (5)*x[1] + (5)*x[2] + (1)*x[3] + (0)*x[4];
s.t. c1: (x[4] <= 2) and (x[3] != 3) ==> (x[1] >= 2);
s.t. c2: abs(x[4] - x[2]) <= 3;
s.t. c3: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c4: (x[4] <= 2) and (x[3] != 3) ==> (x[1] >= 2);
s.t. c5: atleast 2 {i in 1..4} (x[i] >= 2);
