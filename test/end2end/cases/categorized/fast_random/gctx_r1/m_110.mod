var x{1..3} integer >= 0 <= 3;
maximize o: (4)*x[1] + (4)*x[2] + (-4)*x[3];
s.t. c1: (x[1] >= x[1]) and ((x[1] >= x[1]) or not (x[1] >= x[1]));
s.t. c2: (if x[1] >= x[1] then x[1] else 0) + (if x[1] >= x[1] then x[1] else 0) < 1;
s.t. c3: (abs(x[2] - x[2]) <= 0) ==> ((x[1] >= x[1]) or (x[1] >= x[1]));
