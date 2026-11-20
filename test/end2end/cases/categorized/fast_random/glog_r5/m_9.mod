var x{1..3} integer >= -2 <= 2;
maximize o: (-4)*x[1] + (3)*x[2] + (-2)*x[3];
s.t. c1: (x[1] <= -1) <==> (x[3] + x[3] == 2);
s.t. c2: ((x[2] == 1) ==> (x[2] > 1)) and ((x[2] > 1) ==> (x[2] >= 1));
s.t. c3: (if x[2] != 1 then x[1] else x[3]) < 1;
s.t. c4: ((if x[2] != 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
