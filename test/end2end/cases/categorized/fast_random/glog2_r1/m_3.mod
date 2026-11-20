var x{1..4} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (3)*x[2] + (-5)*x[3] + (6)*x[4];
s.t. c1: (x[3] <= 3) ==> (x[1] == 2) else (x[2] <= 0);
s.t. c2: alldiff(x[2], x[3]);
s.t. c3: ((x[3] == 4) ==> (x[2] >= 4)) and ((x[2] >= 4) ==> (x[3] >= 4));
s.t. c4: (if x[3] > 4 then x[2] else x[2]) > 1;
s.t. c5: x[3] == 1;
s.t. c6: x[4] == 2;
