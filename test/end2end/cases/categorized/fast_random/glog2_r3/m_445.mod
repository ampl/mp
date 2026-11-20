var x{1..3} integer >= 0 <= 3;
maximize o: (1)*x[1] + (0)*x[2] + (-5)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: exists {i in 1..3} x[i] == 1;
s.t. c3: (if x[2] < 0 then x[3] else x[1]) <= 1;
s.t. c4: not ((x[3] <= 3) and (x[2] > 2));
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: not ((x[2] != 3) and (x[1] < 1));
