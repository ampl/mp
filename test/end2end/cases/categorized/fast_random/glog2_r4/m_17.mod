var x{1..5} integer >= -2 <= 0;
minimize o: (-1000)*x[1] + (2000)*x[2] + (-999)*x[3] + (-3000)*x[4] + (51)*x[5];
s.t. c1: alldiff(x[2], x[1], x[3]);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) == 3;
s.t. c4: (x[4] < -2) or (x[2] != -2) or (x[3] == 0);
