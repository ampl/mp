var x{1..5} integer >= -3 <= 2;
maximize o: (4)*x[1] + (-5)*x[2] + (-4)*x[3] + (0)*x[4] + (1)*x[5];
s.t. c1: max(x[3], x[2]) - min(x[3], x[4]) >= -3;
s.t. c2: atleast 1 {i in 1..5} (x[i] >= 2);
s.t. c3: (x[5] == -3) <==> (x[3] + x[3] > 2);
