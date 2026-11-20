var x{1..3} integer >= 0 <= 2;

maximize o: -2*x[1] + 0*x[2] - 4*x[3];

s.t. c1: (x[3] <= 0) ==> (x[2] <= 2) else (x[1] > 2);
s.t. c2: atleast 3 {i in 1..3} (x[i] >= 0);
s.t. c3: ((x[3] > 2) ==> (x[1] <= 1)) and ((x[1] <= 1) ==> (x[2] >= 2));
s.t. c4: atmost 2 {i in 1..3} (x[i] >= 0);
