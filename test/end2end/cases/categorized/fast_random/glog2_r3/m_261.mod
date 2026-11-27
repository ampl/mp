var x{1..6} integer >= 0 <= 2;

minimize o: 3*x[1] - 2*x[2] - 2*x[3] + 5*x[4] + x[5] - 3*x[6];

s.t. c1: abs(x[2] - x[5]) < 1;
s.t. c2: atleast 3 {i in 1..6} (x[i] >= -1);
s.t. c3: max(x[5], x[1]) - min(x[3], x[3]) >= -1;
s.t. c4: (x[2] < 0) ==> (x[2] == -1);
