var x{1..3} integer >= -2 <= 2;

minimize o: 4*x[1] - 2*x[2] + 2*x[3];

s.t. c1: atleast 3 {i in 1..3} (x[i] >= 1);
s.t. c2: max(x[3], x[3]) - min(x[1], x[2]) == -2;
s.t. c3: (x[2] != -1) or (x[2] != -1) or (x[1] == -2);
s.t. c4: x[1] == 0;
