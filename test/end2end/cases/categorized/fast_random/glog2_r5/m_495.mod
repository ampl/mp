var x{1..3} integer >= 0 <= 3;

maximize o: 3*x[1] + 0*x[2] - x[3];

s.t. c1: atmost 3 {i in 1..3} (x[i] >= 1);
s.t. c2: exists {i in 1..3} x[i] == 1;
s.t. c3: (x[3] == 2) <==> (x[1] + x[2] < 0);
s.t. c4: abs(x[3] - x[2]) >= 0;
s.t. c5: not ((x[2] != 2) and (x[1] >= 0));
