var x{1..4} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (5)*x[2] + (3)*x[3] + (-2)*x[4];
s.t. c1: (if x[4] > 4 then x[3] else x[1]) != 0;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: (x[1] == 3) and (x[1] != 4) ==> (x[3] >= 3);
s.t. c4: (if x[2] < 3 then x[2] else x[3]) < 4;
