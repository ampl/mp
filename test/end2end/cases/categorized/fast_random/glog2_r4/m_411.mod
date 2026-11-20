var x{1..5} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (-2)*x[2] + (-4)*x[3] + (0)*x[4] + (3)*x[5];
s.t. c1: atleast 1 {i in 1..5} (x[i] >= 3);
s.t. c2: abs(x[1] - x[2]) > 0;
s.t. c3: (x[2] <= 0) and (x[5] < 3) ==> (x[2] >= 0);
