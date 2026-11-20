var x{1..5} integer >= 0 <= 4;
maximize o: (-5)*x[1] + (-1)*x[2] + (1)*x[3] + (-4)*x[4] + (-4)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c3: abs(x[1] - x[5]) > 0;
