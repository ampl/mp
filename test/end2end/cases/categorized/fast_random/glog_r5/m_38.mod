var x{1..3} integer >= 0 <= 2;

minimize o: 3*x[1] + 6*x[2] - 3*x[3];

s.t. c1: atleast 2 {i in 1..3} (x[i] >= 2);
s.t. c2: atleast 3 {i in 1..3} (x[i] >= 1);
s.t. c3: (x[1] > 0) ==> (x[3] != 0);
s.t. c4: not ((x[1] > 0) and (x[1] == 2));
