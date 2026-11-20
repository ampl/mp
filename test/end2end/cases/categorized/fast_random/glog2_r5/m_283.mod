var x{1..5} integer >= -2 <= 2;
maximize o: (-3)*x[1] + (4)*x[2] + (3)*x[3] + (-5)*x[4] + (-2)*x[5];
s.t. c1: abs(x[1] - x[5]) == 4;
s.t. c2: (x[4] >= 0) ==> (x[1] != 2) else (x[3] <= -2);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: ((if x[1] == -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c5: x[5] == -1;
s.t. c6: x[4] == 0;
