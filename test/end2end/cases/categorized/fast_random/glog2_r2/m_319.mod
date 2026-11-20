var x{1..4} integer >= -2 <= 2;

minimize o: 5*x[1] + 5*x[2] - x[3] + 6*x[4];

s.t. c1: (x[3] < 1) <==> (x[3] + x[3] >= -2);
s.t. c2: abs(x[3] - x[2]) < 0;
s.t. c3: abs(x[3] - x[4]) < 1;
s.t. c4: atmost 1 {i in 1..4} (x[i] >= -2);
s.t. c5: (x[3] < 1) <==> (x[3] + x[3] >= -2);
s.t. c6: atmost 1 {i in 1..4} (x[i] >= -2);
