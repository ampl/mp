var x{1..4} integer >= 0 <= 3;
maximize o: (4)*x[1] + (3)*x[2] + (-1)*x[3] + (3)*x[4];
s.t. c1: ((x[1] != 3) ==> (x[2] >= 1)) and ((x[2] >= 1) ==> (x[1] >= 3));
s.t. c2: (x[2] == 0) ==> (x[4] > 2);
s.t. c3: ((if x[1] == 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c4: alldiff(x[2], x[3], x[1], x[4]);
s.t. c5: x[4] == 3;
s.t. c6: x[1] == 0;
