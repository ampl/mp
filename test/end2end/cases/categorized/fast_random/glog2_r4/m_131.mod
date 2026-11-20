var x{1..3} integer >= -2 <= 2;

maximize o: 3*x[1] + 6*x[2] + 2*x[3];

s.t. c1: (x[3] < 3) ==> (x[1] > 3) else (x[2] != 3);
s.t. c2: atmost 1 {i in 1..3} (x[i] >= -3);
s.t. c3: ((x[1] <= 3) ==> (x[3] >= 3)) and ((x[3] >= 3) ==> (x[3] >= 3));
s.t. c4: exists {i in 1..3} x[i] == -3;
