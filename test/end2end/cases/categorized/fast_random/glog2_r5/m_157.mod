var x{1..5} integer >= 0 <= 2;
maximize o: (-3000)*x[1] + (-999)*x[2] + (17)*x[3] + (-2000)*x[4] + (51)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[2] < 0) <==> (x[1] + x[5] <= 2);
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) < 1;
