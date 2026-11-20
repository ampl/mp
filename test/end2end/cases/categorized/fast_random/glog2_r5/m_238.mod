var x{1..4} integer >= 0 <= 1;

minimize o: -x[1] + 6*x[2] + 3*x[3] + x[4];

s.t. c1: atleast 3 {i in 1..4} (x[i] >= 1);
s.t. c2: not ((x[2] <= 1) and (x[3] >= 1));
s.t. c3: (x[3] <= 0) and (x[3] <= 0) ==> (x[3] >= 0);
s.t. c4: ((x[3] > 0) ==> (x[1] == 1)) and ((x[1] == 1) ==> (x[1] >= 0));
