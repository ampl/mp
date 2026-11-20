var x{1..4} integer >= 0 <= 2;
maximize o: (-4)*x[1] + (-3)*x[2] + (1)*x[3] + (5)*x[4];
s.t. c1: atleast 3 {i in 1..4} (x[i] >= 0);
s.t. c2: atleast 1 {i in 1..4} (x[i] >= 0);
s.t. c3: (x[1] >= 2) ==> (x[4] <= 1);
s.t. c4: abs(x[1] - x[3]) >= 1;
s.t. c5: x[2] == 0;
s.t. c6: x[3] == 0;
