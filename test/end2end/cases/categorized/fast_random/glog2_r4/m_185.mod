var x{1..3} integer >= 0 <= 4;
maximize o: (0)*x[1] + (-4)*x[2] + (2)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 4;
s.t. c2: atleast 3 {i in 1..3} (x[i] >= 0);
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: max(x[2], x[3]) - min(x[2], x[1]) < 3;
s.t. c5: (x[3] < 3) ==> (x[2] <= 0);
s.t. c6: numberof 3 in (x[1], x[2], x[3]) >= 2;
s.t. c7: (x[2] == 3) ==> (x[1] > 2);
