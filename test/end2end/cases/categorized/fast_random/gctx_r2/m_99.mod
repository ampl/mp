var x{1..3} integer >= -1 <= 2;
maximize o: (-1)*x[1] + (-3)*x[2] + (-1)*x[3];
s.t. c1: (abs(x[2] - x[1]) > 2) <==> (min(x[3], x[1]) >= -1);
s.t. c2: (abs(x[2] - x[1]) > 2) or (min(x[3], x[1]) >= -1);
s.t. c3: (abs(x[2] - x[1]) > 2) or (min(x[3], x[1]) >= -1);
s.t. c4: (abs(x[2] - x[1]) > 2) <==> (abs(x[2] - x[1]) > 2);
