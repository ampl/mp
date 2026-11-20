var x{1..6} integer >= 0 <= 2;

minimize o: -2*x[1] + 3*x[2] - 4*x[3] + 3*x[4] + 3*x[5] + 5*x[6];

s.t. c1: atleast 1 {i in 1..6} (x[i] >= 2);
s.t. c2: ((x[1] >= 0) ==> (x[2] <= 1)) and ((x[2] <= 1) ==> (x[4] >= 0));
s.t. c3: max(x[2], x[3]) - min(x[6], x[2]) > 0;
s.t. c4: atleast 1 {i in 1..6} (x[i] >= 2);
s.t. c5: max(x[2], x[3]) - min(x[6], x[2]) > 0;
