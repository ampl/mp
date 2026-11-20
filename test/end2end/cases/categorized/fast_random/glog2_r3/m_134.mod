var x{1..3} integer >= -2 <= 2;
maximize o: (6)*x[1] + (-5)*x[2] + (-3)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: atleast 2 {i in 1..3} (x[i] >= 2);
