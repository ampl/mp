var x{1..5} integer >= 0 <= 2;

minimize o: x[1] + 3*x[2] - 2*x[3] - 2*x[4] + 0*x[5];

s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (x[2] < 1) ==> (x[5] != 1);
s.t. c3: atleast 2 {i in 1..5} (x[i] >= 1);
s.t. c4: (x[1] > 1) ==> (x[3] >= 2) else (x[4] < 1);
s.t. c5: abs(x[1] - x[3]) > 0;
s.t. c6: max(x[4], x[2]) - min(x[4], x[3]) <= 1;
s.t. c7: abs(x[4] - x[4]) >= 1;
