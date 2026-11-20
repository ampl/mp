var x{1..4} integer >= -2 <= 0;
maximize o: (-1)*x[1] + (1)*x[2] + (6)*x[3] + (-1)*x[4];
s.t. c1: ((x[1] != -2) ==> (x[4] < -1)) and ((x[4] < -1) ==> (x[2] >= -2));
s.t. c2: ((if x[1] <= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: numberof -1 in (x[1], x[2], x[3], x[4]) > 0;
s.t. c4: exists {i in 1..4} x[i] == -1;
s.t. c5: exists {i in 1..4} x[i] == 0;
