var x{1..3} integer >= -2 <= 0;
maximize o: (-2000)*x[1] + (17)*x[2] + (-999)*x[3];
s.t. c1: atmost 2 {i in 1..3} (x[i] >= 0);
