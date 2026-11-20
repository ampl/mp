var x{1..3} integer >= -3 <= 2;

minimize o: 0*x[1] - 2*x[2] + 3*x[3];

s.t. c1: ((x[3] <= 0) ==> (x[3] > 1)) and ((x[3] > 1) ==> (x[2] >= 0));
s.t. c2: atleast 3 {i in 1..3} (x[i] >= -2);
s.t. c3: not ((x[3] > 1) and (x[2] >= 0));
