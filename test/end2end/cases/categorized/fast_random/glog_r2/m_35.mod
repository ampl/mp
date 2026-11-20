var x{1..6} integer >= -2 <= 1;
minimize o: (-1)*x[1] + (3)*x[2] + (-5)*x[3] + (-2)*x[4] + (1)*x[5] + (3)*x[6];
s.t. c1: not ((x[5] == 1) and (x[1] < -2));
s.t. c2: (x[6] < -1) ==> (x[2] >= 1) else (x[6] > -2);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
