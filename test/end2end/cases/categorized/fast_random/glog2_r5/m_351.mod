var x{1..3} integer >= 0 <= 2;
minimize o: (5)*x[1] + (-4)*x[2] + (1)*x[3];
s.t. c1: numberof 0 in (x[1], x[2], x[3]) > 1;
s.t. c2: ((x[3] > 2) ==> (x[3] == 0)) and ((x[3] == 0) ==> (x[3] >= 2));
s.t. c3: exists {i in 1..3} x[i] == 0;
