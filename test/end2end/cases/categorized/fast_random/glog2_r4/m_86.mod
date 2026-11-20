var x{1..3} integer >= -2 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: (x[1] > -1) or (x[3] < 1) or (x[3] > 1);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: numberof -2 in (x[1], x[2], x[3]) < 3;
s.t. c4: ((x[3] >= -1) ==> (x[1] > 0)) and ((x[1] > 0) ==> (x[3] >= -1));
