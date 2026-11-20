var x{1..3} integer >= 0 <= 2;

maximize o: -3*x[1] + x[2] - 5*x[3];

s.t. c1: atmost 3 {i in 1..3} (x[i] >= 2);
s.t. c2: (if x[2] >= 1 then x[3] else x[3]) <= 0;
s.t. c3: (x[3] == 0) or (x[2] < 2) or (x[1] <= 0);
s.t. c4: (if x[3] > 1 then x[2] else x[2]) == 1;
