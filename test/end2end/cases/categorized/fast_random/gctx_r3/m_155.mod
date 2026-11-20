var x{1..4} integer >= 0 <= 3;
maximize o: (4)*x[1] + (4)*x[2] + (2)*x[3] + (0)*x[4] + (2)*(if max(x[2], x[2]) >= 2 then 1 else 0);
s.t. c1: (abs(x[3] - x[4]) >= 0) ==> ((abs(x[3] - x[4]) >= 0) or (x[1] < x[3]));
s.t. c2: (if abs(x[3] - x[4]) >= 0 then 1 else 0) + (if max(x[2], x[2]) >= 2 then 1 else 0) + (if max(x[2], x[2]) >= 2 then 1 else 0) >= 2;
s.t. c3: (abs(x[3] - x[4]) >= 0) <==> (abs(x[3] - x[4]) >= 0);
