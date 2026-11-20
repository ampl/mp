var x{1..5} integer >= -2 <= 0;

maximize o: 5*x[1] + 3*x[2] - x[3] - x[4] + 4*x[5];

s.t. c1: abs(x[1] - x[4]) != 0;
s.t. c2: exists {i in 1..5} x[i] == 0;
s.t. c3: (if x[2] != -2 then x[2] else x[5]) >= 0;
s.t. c4: x[2] == -1;
s.t. c5: x[1] == 0;
