var x{1..5} integer >= -2 <= 1;
minimize o: (4)*x[1] + (1)*x[2] + (2)*x[3] + (4)*x[4] + (-5)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (if x[5] > -1 then x[5] else x[5]) >= 0;
s.t. c3: not ((x[3] != -1) and (x[2] != 0));
s.t. c4: x[5] == -2;
s.t. c5: x[4] == -2;
