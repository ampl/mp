var x{1..6} integer >= -2 <= 1;
minimize o: (4)*x[1] + (-5)*x[2] + (2)*x[3] + (-5)*x[4] + (6)*x[5] + (-3)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[5] <= -2) or (x[3] > -1) or (x[5] < 1);
s.t. c3: x[2] == 1;
