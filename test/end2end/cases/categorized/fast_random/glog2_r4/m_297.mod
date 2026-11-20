var x{1..5} integer >= 0 <= 2;
minimize o: (-666)*x[1] + (51)*x[2] + (-3000)*x[3] + (1000)*x[4] + (17)*x[5];
s.t. c1: atmost 3 {i in 1..5} (x[i] >= 1);
