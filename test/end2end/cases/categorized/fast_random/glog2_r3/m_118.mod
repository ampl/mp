var x{1..7} integer >= 0 <= 1;

minimize o: 2*x[1] - 4*x[2] - 2*x[3] - 2*x[4] + 4*x[5] + 0*x[6] - 3*x[7];

s.t. c1: ((x[3] != 1) ==> (x[1] == 0)) and ((x[1] == 0) ==> (x[7] >= 1));
s.t. c2: exists {i in 1..7} x[i] == 1;
s.t. c3: atleast 1 {i in 1..7} (x[i] >= 1);
