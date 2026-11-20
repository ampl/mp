var x{1..5} integer >= -1 <= 1;
maximize o: (5)*x[1] + (5)*x[2] + (-5)*x[3] + (1)*x[4] + (-3)*x[5] + (-1)*(if max(x[1], x[1]) > 1 then 1 else 0);
s.t. c1: (abs(x[5] - x[5]) >= 1) <==> (max(x[1], x[1]) > 1);
s.t. c2: (max(x[1], x[1]) > 1) and ((max(x[1], x[1]) > 1) or not (x[3] != -1));
