var x{1..3} integer >= 0 <= 3;
maximize o: (2)*x[1] + (1)*x[2] + (1)*x[3];
s.t. c1: numberof 1 in (x[1], x[2], x[3]) < 1;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: ((if x[3] != 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c4: (if x[3] >= 0 then x[3] else x[1]) > 0;
