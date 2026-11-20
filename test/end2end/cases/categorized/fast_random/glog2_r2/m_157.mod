var x{1..5} integer >= 0 <= 4;
minimize o: (-2000)*x[1] + (34)*x[2] + (-999)*x[3] + (-666)*x[4] + (1000)*x[5];
s.t. c1: atmost 1 {i in 1..5} (x[i] >= 2);
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: max(x[1], x[5]) - min(x[4], x[2]) > 0;
