var x{1..4} integer >= 0 <= 3;

minimize o: -x[1] + 0*x[2] - 2*x[3] - 5*x[4];

s.t. c1: (x[3] > 2) and (x[4] > 2) ==> (x[4] >= 2);
s.t. c2: atmost 3 {i in 1..4} (x[i] >= 3);
s.t. c3: atleast 3 {i in 1..4} (x[i] >= 2);
