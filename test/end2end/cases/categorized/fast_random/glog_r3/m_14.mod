var x{1..3} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (-1)*x[2] + (1)*x[3];
s.t. c1: not ((x[1] != 1) and (x[2] < 3));
s.t. c2: (x[3] >= 2) ==> (x[2] >= 1) else (x[1] > 0);
s.t. c3: atleast 3 {i in 1..3} (x[i] >= 2);
