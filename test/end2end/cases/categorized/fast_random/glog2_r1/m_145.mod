var x{1..5} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (5)*x[2] + (-3)*x[3] + (-2)*x[4] + (-2)*x[5];
s.t. c1: (if x[1] < 3 then x[3] else x[5]) != 4;
s.t. c2: (x[2] >= 0) ==> (x[2] != 4) else (x[3] == 2);
s.t. c3: alldiff(x[1], x[4], x[3], x[5]);
s.t. c4: (if x[2] < 1 then x[4] else x[5]) <= 4;
s.t. c5: max(x[3], x[3]) - min(x[5], x[3]) <= 1;
s.t. c6: (x[1] <= 0) <==> (x[4] + x[5] <= 3);
s.t. c7: (x[4] == 0) ==> (x[4] >= 4);
