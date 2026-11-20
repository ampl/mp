var x{1..5} integer >= -2 <= 2;
minimize o: (-1)*x[1] + (1)*x[2] + (0)*x[3] + (-4)*x[4] + (6)*x[5];
s.t. c1: (x[4] != 2) <==> (x[1] + x[1] < 2);
s.t. c2: (x[1] == 0) ==> (x[4] != -2);
s.t. c3: atmost 3 {i in 1..5} (x[i] >= 1);
s.t. c4: atleast 1 {i in 1..5} (x[i] >= -2);
s.t. c5: (x[4] != 2) <==> (x[1] + x[1] < 2);
s.t. c6: atleast 1 {i in 1..5} (x[i] >= -2);
