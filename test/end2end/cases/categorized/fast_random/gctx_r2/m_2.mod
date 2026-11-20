var x{1..4} integer >= 0 <= 2;
maximize o: (4)*x[1] + (5)*x[2] + (2)*x[3] + (-4)*x[4] + (0)*(if x[2] + x[1] < 0 then 1 else 0);
s.t. c1: (x[2] + x[1] < 0) and ((x[2] + x[1] < 0) or not (abs(x[1] - x[3]) == 0));
s.t. c2: (abs(x[4] - x[2]) >= 0) and ((abs(x[4] - x[2]) >= 0) or not (abs(x[4] - x[2]) >= 0));
s.t. c3: (x[2] + x[1] < 0) ==> ((abs(x[4] - x[2]) >= 0) or (x[2] + x[1] < 0));
