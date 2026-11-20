var x{1..6} integer >= 0 <= 3;
minimize o: (-999)*x[1] + (-1000)*x[2] + (17)*x[3] + (2000)*x[4] + (17)*x[5] + (17)*x[6];
s.t. c1: atleast 1 {i in 1..6} (x[i] >= 1);
