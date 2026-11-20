var x{1..3} integer >= 0 <= 3;

maximize o: x[1] + 5*x[2] - 3*x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= 3);
s.t. c2: (x[1] >= 1) or (x[1] <= 2) or (x[2] < 3);
s.t. c3: (x[3] == 3) and (x[1] <= 2) ==> (x[3] >= 3);
