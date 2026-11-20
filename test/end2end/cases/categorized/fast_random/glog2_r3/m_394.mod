var x{1..3} integer >= -2 <= 3;
minimize o: (-4)*x[1] + (-3)*x[2] + (-5)*x[3];
s.t. c1: (x[2] < 1) <==> (x[3] + x[2] >= -2);
s.t. c2: (if x[2] < 3 then x[3] else x[2]) < 1;
s.t. c3: (if x[2] == 3 then x[2] else x[1]) < -1;
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
