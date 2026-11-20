var x{1..6} integer >= 0 <= 3;

maximize o: -2*x[1] + 2*x[2] - 4*x[3] - 2*x[4] + 0*x[5] - 2*x[6];

s.t. c1: not ((x[1] != -1) and (x[3] <= 0));
s.t. c2: (x[4] > 4) and (x[1] == -1) ==> (x[4] >= 4);
s.t. c3: atleast 3 {i in 1..6} (x[i] >= 3);
s.t. c4: abs(x[3] - x[2]) <= 1;
