var x{1..3} integer >= -2 <= 1;
minimize o: (-3)*x[1] + (-4)*x[2] + (-1)*x[3];
s.t. c1: exists {i in 1..3} x[i] == -1;
s.t. c2: abs(x[1] - x[2]) < 3;
s.t. c3: (x[3] != -1) <==> (x[1] + x[3] != -2);
s.t. c4: (x[2] > 1) and (x[1] <= -1) ==> (x[1] >= 1);
