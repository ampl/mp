var x{1..3} integer >= 0 <= 2;
maximize o: (1)*x[1] + (5)*x[2] + (-2)*x[3];
s.t. c1: (max(x[3], x[2]) <= 0) or (x[3] + x[2] <= 2);
s.t. c2: (max(x[3], x[2]) <= 0) <==> (x[3] < 0);
s.t. c3: (x[3] < 0) and ((x[3] + x[2] <= 2) or not (x[3] + x[2] <= 2));
