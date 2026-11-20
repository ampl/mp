var x{1..5} integer >= 0 <= 2;
maximize o: (4)*x[1] + (1)*x[2] + (-1)*x[3] + (6)*x[4] + (-3)*x[5] + (-5)*(if max(x[1], x[2]) != 0 then 1 else 0);
s.t. c1: (abs(x[2] - x[5]) < 0) or (max(x[1], x[2]) != 1);
s.t. c2: (max(x[1], x[2]) != 1) ==> (max(x[1], x[2]) != 2);
