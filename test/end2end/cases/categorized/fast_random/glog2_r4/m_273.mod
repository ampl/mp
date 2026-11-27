var x{1..6} integer >= -2 <= 0;

minimize o: -4*x[1] - 4*x[2] - 5*x[3] + 0*x[4] - 5*x[5] + 6*x[6];

s.t. c1: exists {i in 1..6} x[i] == -2;
s.t. c2: (x[6] <= -2) ==> (x[3] != -1);
s.t. c3: max(x[5], x[5]) - min(x[3], x[2]) == -1;
s.t. c4: atleast 3 {i in 1..6} (x[i] >= -2);
s.t. c5: x[4] == -2;
s.t. c6: x[6] == -1;
