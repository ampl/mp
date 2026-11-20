var x{1..3} integer >= 0 <= 2;
minimize o: (-2)*x[1] + (-1)*x[2] + (-4)*x[3];
s.t. c1: (x[2] != 0) and (x[1] >= 0) ==> (x[1] >= 0);
s.t. c2: (x[3] >= 1) ==> (x[2] < 0);
s.t. c3: (if x[3] >= 2 then x[3] else x[3]) <= 2;
s.t. c4: alldiff(x[3], x[1], x[2]);
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: ((x[1] == 1) ==> (x[3] <= 1)) and ((x[3] <= 1) ==> (x[3] >= 1));
s.t. c7: numberof 1 in (x[1], x[2], x[3]) < 3;
