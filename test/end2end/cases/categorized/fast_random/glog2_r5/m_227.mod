var x{1..5} integer >= -2 <= 2;
maximize o: (34)*x[1] + (51)*x[2] + (17)*x[3] + (51)*x[4] + (-999)*x[5];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) > 3;
s.t. c2: exists {i in 1..5} x[i] == 2;
