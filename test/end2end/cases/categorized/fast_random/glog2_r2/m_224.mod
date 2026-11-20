var x{1..5} integer >= -2 <= 3;

maximize o: 3*x[1] + 4*x[2] + 0*x[3] + 6*x[4] + 2*x[5];

s.t. c1: (if x[3] < 3 then x[2] else x[5]) < 3;
s.t. c2: atmost 1 {i in 1..5} (x[i] >= 3);
