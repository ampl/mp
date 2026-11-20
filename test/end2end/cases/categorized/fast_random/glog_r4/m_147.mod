var x{1..5} integer >= -2 <= 1;
maximize o: (3)*x[1] + (1)*x[2] + (4)*x[3] + (4)*x[4] + (1)*x[5];
s.t. c1: max(x[2], x[5]) - min(x[3], x[3]) != 0;
s.t. c2: (x[4] > 1) <==> (x[5] + x[5] < 1);
s.t. c3: not ((x[1] < -1) and (x[1] != 1));
s.t. c4: atleast 2 {i in 1..5} (x[i] >= 0);
