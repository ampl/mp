var x{1..6} integer >= -2 <= 2;

minimize o: 6*x[1] + x[2] + x[3] + 6*x[4] - 5*x[5] + 2*x[6];

s.t. c1: abs(x[1] - x[3]) != 1;
s.t. c2: (x[3] > 0) or (x[5] == 2) or (x[4] < 2);
s.t. c3: ((x[1] >= 2) ==> (x[6] <= -2)) and ((x[6] <= -2) ==> (x[2] >= 2));
s.t. c4: exists {i in 1..6} x[i] == 0;
