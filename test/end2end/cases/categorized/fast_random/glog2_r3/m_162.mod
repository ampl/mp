var x{1..3} integer >= -2 <= 0;

minimize o: 5*x[1] + 6*x[2] + x[3];

s.t. c1: exists {i in 1..3} x[i] == -1;
s.t. c2: abs(x[3] - x[3]) < 2;
s.t. c3: (x[1] <= -1) <==> (x[3] + x[1] == -1);
s.t. c4: (x[2] < 0) ==> (x[2] > 0) else (x[2] == 0);
