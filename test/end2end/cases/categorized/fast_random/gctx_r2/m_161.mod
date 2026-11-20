var x{1..5} integer >= -1 <= 1;
minimize o: (4)*x[1] + (5)*x[2] + (4)*x[3] + (-4)*x[4] + (4)*x[5];
s.t. c1: not ((max(x[4], x[1]) >= -1) and (abs(x[2] - x[2]) == 2));
s.t. c2: not (abs(x[2] - x[2]) == 2);
s.t. c3: (abs(x[2] - x[2]) == 2) <==> (max(x[4], x[1]) >= -1);
s.t. c4: (abs(x[2] - x[2]) == 2) <==> (max(x[4], x[1]) >= -1);
