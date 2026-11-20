var x{1..3} integer >= 0 <= 3;

maximize o: 2000*x[1] + 17*x[2] + 51*x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= 3);
s.t. c2: max(x[3], x[2]) - min(x[2], x[2]) < 2;
s.t. c3: not ((x[2] < 3) and (x[3] > 1));
s.t. c4: max(x[1], x[1]) - min(x[1], x[1]) < 3;
