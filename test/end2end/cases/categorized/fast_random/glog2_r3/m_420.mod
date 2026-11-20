var x{1..5} integer >= -2 <= 2;
minimize o: (-3)*x[1] + (-1)*x[2] + (6)*x[3] + (2)*x[4] + (-4)*x[5];
s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) < 3;
s.t. c2: atleast 3 {i in 1..5} (x[i] >= -2);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: (x[2] < -1) or (x[3] != 1) or (x[2] == -2);
