var x{1..6} integer >= -2 <= 2;
minimize o: (6)*x[1] + (-4)*x[2] + (2)*x[3] + (-1)*x[4] + (-5)*x[5] + (-1)*x[6];
s.t. c1: (x[3] > -2) <==> (x[1] + x[5] != 1);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: (if x[4] == 0 then x[6] else x[3]) > -2;
s.t. c4: (x[3] > -2) <==> (x[1] + x[5] != 1);
s.t. c5: (if x[4] == 0 then x[6] else x[3]) > -2;
