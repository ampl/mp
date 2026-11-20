var x{1..5} integer >= 0 <= 4;
maximize o: (-2)*x[1] + (4)*x[2] + (3)*x[3] + (3)*x[4] + (3)*x[5];
s.t. c1: not ((x[2] < 1) and (x[2] != 4));
s.t. c2: ((x[4] > 4) ==> (x[2] == 1)) and ((x[2] == 1) ==> (x[3] >= 4));
s.t. c3: (x[3] > 4) ==> (x[3] == 3) else (x[4] == 2);
s.t. c4: exists {i in 1..5} x[i] == 0;
s.t. c5: exists {i in 1..5} x[i] == 3;
s.t. c6: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) < 3;
s.t. c7: max(x[1], x[3]) - min(x[1], x[2]) > 2;
