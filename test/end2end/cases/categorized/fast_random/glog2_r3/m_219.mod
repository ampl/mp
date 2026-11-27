var x{1..3} integer >= 0 <= 3;

minimize o: -x[1] + 4*x[2] + 2*x[3];

s.t. c1: (x[2] == 3) <==> (x[3] + x[3] < 0);
s.t. c2: exists {i in 1..3} x[i] == 2;
s.t. c3: max(x[2], x[1]) - min(x[3], x[2]) == 0;
s.t. c4: (x[2] == 3) <==> (x[3] + x[3] < 0);
s.t. c5: max(x[2], x[1]) - min(x[3], x[2]) == 0;
