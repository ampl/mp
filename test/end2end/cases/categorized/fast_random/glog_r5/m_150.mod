var x{1..5} integer >= 0 <= 2;

minimize o: -2*x[1] - 3*x[2] + 2*x[3] + 2*x[4] - 4*x[5];

s.t. c1: (x[1] == 0) <==> (x[2] + x[2] < 1);
s.t. c2: atleast 3 {i in 1..5} (x[i] >= 0);
s.t. c3: atmost 2 {i in 1..5} (x[i] >= 2);
s.t. c4: (x[3] == 0) ==> (x[3] != 0) else (x[4] != 0);
