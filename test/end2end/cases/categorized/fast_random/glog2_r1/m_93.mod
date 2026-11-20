var x{1..5} integer >= -2 <= 1;
minimize o: (-4)*x[1] + (6)*x[2] + (4)*x[3] + (3)*x[4] + (0)*x[5];
s.t. c1: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c2: atleast 3 {i in 1..5} (x[i] >= 1);
s.t. c3: ((x[1] > 1) ==> (x[4] >= -2)) and ((x[4] >= -2) ==> (x[5] >= 1));
s.t. c4: x[1] == -2;
