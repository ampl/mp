var x{1..3} integer >= 0 <= 3;
minimize o: (4)*x[1] + (5)*x[2] + (0)*x[3];
s.t. c1: atmost 1 {i in 1..3} (x[i] >= 2);
s.t. c2: (x[3] < 2) and (x[1] == 1) ==> (x[2] >= 2);
s.t. c3: not ((x[3] == 3) and (x[2] < 0));
s.t. c4: ((x[3] > 3) ==> (x[1] > 2)) and ((x[1] > 2) ==> (x[1] >= 3));
s.t. c5: abs(x[2] - x[2]) != 0;
s.t. c6: forall {i in 1..2} x[i] != x[i+1];
