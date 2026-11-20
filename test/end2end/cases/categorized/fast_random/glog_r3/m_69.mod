var x{1..3} integer >= -2 <= 2;

maximize o: -3*x[1] - 4*x[2] + 3*x[3];

s.t. c1: (if x[2] <= 0 then x[3] else x[2]) > 0;
s.t. c2: atmost 2 {i in 1..3} (x[i] >= -1);
