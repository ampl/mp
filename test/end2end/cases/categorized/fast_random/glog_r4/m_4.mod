var x{1..4} integer >= 0 <= 3;
minimize o: (3)*x[1] + (3)*x[2] + (-1)*x[3] + (-4)*x[4];
s.t. c1: (x[1] != 3) and (x[4] != 2) ==> (x[3] >= 3);
s.t. c2: atleast 1 {i in 1..4} (x[i] >= 2);
s.t. c3: abs(x[4] - x[4]) <= 0;
