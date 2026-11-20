var x{1..5} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (-1)*x[2] + (4)*x[3] + (-2)*x[4] + (4)*x[5] + (-2)*(if abs(x[1] - x[5]) < 1 then 1 else 0);
s.t. c1: (x[4] == 2) ==> (x[5] <= 1) else (x[4] == 2);
s.t. c2: (abs(x[1] - x[5]) < 1) or (x[5] <= 1);
