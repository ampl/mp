var x{1..3} integer >= 0 <= 2;
minimize o: (3)*x[1] + (1)*x[2] + (2)*x[3];
s.t. c1: (x[1] == 3) ==> (x[1] > 0);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (if x[1] <= -1 then x[1] else x[2]) != 0;
s.t. c4: ((x[1] != -1) ==> (x[2] < 2)) and ((x[2] < 2) ==> (x[2] >= -1));
