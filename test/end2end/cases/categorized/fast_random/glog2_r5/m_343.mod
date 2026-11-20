var x{1..4} integer >= 0 <= 3;
minimize o: (2)*x[1] + (2)*x[2] + (-5)*x[3] + (0)*x[4];
s.t. c1: (x[1] >= 0) <==> (x[2] + x[1] < 1);
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4]) < 2;
s.t. c3: atleast 1 {i in 1..4} (x[i] >= 1);
s.t. c4: x[3] == 2;
s.t. c5: x[1] == 2;
