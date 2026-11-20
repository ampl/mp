var x{1..3} integer >= 0 <= 3;
maximize o: (3)*x[1] + (3)*x[2] + (4)*x[3] + (-3)*(if abs(x[1] - x[2]) <= 3 then 1 else 0);
s.t. c1: (abs(x[1] - x[2]) <= 3) <==> (abs(x[1] - x[2]) <= 3);
s.t. c2: (abs(x[1] - x[2]) <= 3) ==> (abs(x[1] - x[2]) <= 3);
s.t. c3: (if abs(x[1] - x[2]) <= 3 then x[1] else 0) + (if abs(x[1] - x[2]) <= 3 then x[1] else 0) < 2;
