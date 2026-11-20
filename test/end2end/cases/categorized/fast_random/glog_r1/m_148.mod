var x{1..3} integer >= 0 <= 4;

maximize o: 2*x[1] + 6*x[2] + 4*x[3];

s.t. c1: atleast 2 {i in 1..3} (x[i] >= 2);
s.t. c2: atmost 3 {i in 1..3} (x[i] >= 0);
s.t. c3: (x[2] <= 0) or (x[1] != 2) or (x[3] == 0);
s.t. c4: (x[1] != 3) <==> (x[2] + x[1] != 3);
