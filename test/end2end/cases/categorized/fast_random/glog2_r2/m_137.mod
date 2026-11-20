var x{1..5} integer >= -2 <= 0;

minimize o: 3000*x[1] + 34*x[2] - 2000*x[3] + 34*x[4] + 34*x[5];

s.t. c1: max(x[3], x[3]) - min(x[4], x[3]) > -2;
s.t. c2: atmost 2 {i in 1..5} (x[i] >= -1);
s.t. c3: (x[2] <= -2) <==> (x[2] + x[5] != -2);
s.t. c4: max(x[3], x[1]) - min(x[2], x[4]) <= -2;
