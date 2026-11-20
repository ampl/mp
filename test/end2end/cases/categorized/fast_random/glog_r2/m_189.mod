var x{1..5} integer >= 0 <= 4;
minimize o: (4)*x[1] + (3)*x[2] + (-4)*x[3] + (3)*x[4] + (6)*x[5];
s.t. c1: (x[5] == 4) or (x[5] == 4) or (x[3] < 3);
s.t. c2: not ((x[2] != 3) and (x[1] > 0));
s.t. c3: alldiff(x[1], x[4], x[5], x[3]);
s.t. c4: atmost 2 {i in 1..5} (x[i] >= 1);
