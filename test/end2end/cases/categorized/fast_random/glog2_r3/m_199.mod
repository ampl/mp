var x{1..5} integer >= -2 <= 1;

minimize o: 6*x[1] - 4*x[2] + 0*x[3] - 3*x[4] - 3*x[5];

s.t. c1: (x[5] <= -1) or (x[5] < -1) or (x[5] <= -2);
s.t. c2: (x[5] <= -1) ==> (x[3] == -2);
s.t. c3: atleast 1 {i in 1..5} (x[i] >= 0);
s.t. c4: abs(x[4] - x[3]) > 0;
s.t. c5: (x[5] <= -1) or (x[5] < -1) or (x[5] <= -2);
s.t. c6: abs(x[4] - x[3]) > 0;
