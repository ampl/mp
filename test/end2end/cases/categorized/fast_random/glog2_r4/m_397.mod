var x{1..4} integer >= -2 <= 0;

maximize o: 2000*x[1] + 3000*x[2] + 34*x[3] - 999*x[4];

s.t. c1: (if x[1] != 0 then x[1] else x[1]) > -1;
s.t. c2: exists {i in 1..4} x[i] == -1;
