var x{1..3} integer >= 0 <= 2;
minimize o: (0)*x[1] + (6)*x[2] + (0)*x[3];
s.t. c1: atleast 1 {i in 1..3} (x[i] >= 2);
s.t. c2: (x[2] > 2) ==> (x[3] <= 2);
s.t. c3: not ((x[1] < 1) and (x[3] != 1));
s.t. c4: numberof 0 in (x[1], x[2], x[3]) < 3;
