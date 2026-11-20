var x{1..6} integer >= 0 <= 3;
maximize o: (-333)*x[1] + (51)*x[2] + (17)*x[3] + (17)*x[4] + (17)*x[5] + (2000)*x[6];
s.t. c1: atleast 3 {i in 1..6} (x[i] >= 3);
s.t. c2: alldiff(x[1], x[5], x[6]);
