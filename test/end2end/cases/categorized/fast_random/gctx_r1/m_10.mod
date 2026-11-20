var x{1..5} integer >= -1 <= 1;
maximize o: (4)*x[1] + (-5)*x[2] + (2)*x[3] + (3)*x[4] + (2)*x[5];
s.t. c1: (x[3] >= x[1]) ==> (max(x[5], x[3]) != 1);
s.t. c2: (if x[5] <= x[3] then x[5] else 0) + (if x[3] >= x[1] then x[1] else 0) > 0;
s.t. c3: (x[5] <= x[3]) and ((max(x[5], x[3]) != 1) or not (x[5] <= x[3]));
s.t. c4: (max(x[5], x[3]) != 1) or (max(x[5], x[3]) != 1);
