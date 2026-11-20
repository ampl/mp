var x{1..3} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (3)*x[2] + (4)*x[3];
s.t. c1: (if x[3] < 1 then x[3] else x[2]) != 0;
s.t. c2: atmost 3 {i in 1..3} (x[i] >= 0);
s.t. c3: numberof 4 in (x[1], x[2], x[3]) != 2;
s.t. c4: (x[2] >= 3) and (x[3] != 1) ==> (x[3] >= 3);
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: forall {i in 1..2} x[i] != x[i+1];
s.t. c7: (x[2] < 1) and (x[2] == 0) ==> (x[3] >= 1);
