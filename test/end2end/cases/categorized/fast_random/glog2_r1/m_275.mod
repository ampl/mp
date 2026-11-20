var x{1..3} integer >= 0 <= 4;

maximize o: 2*x[1] + 5*x[2] + 3*x[3];

s.t. c1: not ((x[1] < 4) and (x[2] <= 4));
s.t. c2: (x[1] > 4) or (x[2] <= 0) or (x[3] != 1);
s.t. c3: not ((x[1] < 2) and (x[3] == 0));
s.t. c4: atmost 1 {i in 1..3} (x[i] >= 3);
s.t. c5: abs(x[1] - x[2]) <= 3;
