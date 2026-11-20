var x{1..3} integer >= 0 <= 2;
minimize o: (4)*x[1] + (2)*x[2] + (1)*x[3];
s.t. c1: atleast 3 {i in 1..3} (x[i] >= 2);
s.t. c2: (x[3] != 0) or (x[2] >= 1) or (x[3] > 0);
s.t. c3: atleast 3 {i in 1..3} (x[i] >= 2);
s.t. c4: (x[3] != 0) or (x[2] >= 1) or (x[3] > 0);
