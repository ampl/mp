var x{1..3} integer >= 0 <= 2;

minimize o: 0*x[1] + 6*x[2] - 4*x[3];

s.t. c1: (x[3] != 1) ==> (x[2] > 2) else (x[2] > 0);
s.t. c2: atleast 1 {i in 1..3} (x[i] >= 0);
s.t. c3: (x[3] != 1) ==> (x[2] > 2) else (x[2] > 0);
s.t. c4: atleast 1 {i in 1..3} (x[i] >= 0);
