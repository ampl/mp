var x{1..3} integer >= 0 <= 4;
minimize o: (-3)*x[1] + (-5)*x[2] + (3)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: alldiff(x[2], x[3]);
s.t. c3: ((x[2] == 4) ==> (x[3] >= 3)) and ((x[3] >= 3) ==> (x[1] >= 4));
s.t. c4: numberof 2 in (x[1], x[2], x[3]) <= 0;
