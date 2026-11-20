var x{1..3} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (-1)*x[2] + (4)*x[3];
s.t. c1: not ((x[1] == 3) and (x[2] != 1));
s.t. c2: exists {i in 1..3} x[i] == 0;
s.t. c3: atleast 2 {i in 1..3} (x[i] >= 0);
s.t. c4: not ((x[1] == 3) and (x[2] != 1));
s.t. c5: atleast 2 {i in 1..3} (x[i] >= 0);
