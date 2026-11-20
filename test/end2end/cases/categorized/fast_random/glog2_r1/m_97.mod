var x{1..5} integer >= -2 <= 2;
minimize o: (2000)*x[1] + (-1000)*x[2] + (-1000)*x[3] + (-333)*x[4] + (17)*x[5];
s.t. c1: not ((x[1] >= 2) and (x[1] > -2));
s.t. c2: (x[1] == 1) ==> (x[2] > 0) else (x[3] < 1);
s.t. c3: atleast 2 {i in 1..5} (x[i] >= -2);
s.t. c4: atmost 2 {i in 1..5} (x[i] >= 0);
