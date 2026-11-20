var x{1..4} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (6)*x[2] + (3)*x[3] + (2)*x[4];
s.t. c1: atleast 2 {i in 1..4} (x[i] >= 1);
s.t. c2: exists {i in 1..4} x[i] == 2;
s.t. c3: x[2] == 0;
