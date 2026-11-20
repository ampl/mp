var x{1..3} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (-1)*x[2] + (-1)*x[3];
s.t. c1: ((if x[2] == 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c2: ((x[3] >= 1) ==> (x[1] <= 3)) and ((x[1] <= 3) ==> (x[2] >= 1));
s.t. c3: numberof 4 in (x[1], x[2], x[3]) == 1;
s.t. c4: (x[3] > 0) <==> (x[2] + x[1] < 0);
s.t. c5: numberof 3 in (x[1], x[2], x[3]) != 2;
s.t. c6: forall {i in 1..2} x[i] != x[i+1];
s.t. c7: (x[1] <= 1) <==> (x[1] + x[2] != 2);
