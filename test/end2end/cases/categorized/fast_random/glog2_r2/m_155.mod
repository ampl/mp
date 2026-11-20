var x{1..4} integer >= -2 <= 1;
maximize o: (-2)*x[1] + (-1)*x[2] + (1)*x[3] + (-2)*x[4];
s.t. c1: (x[2] < -1) and (x[4] != -1) ==> (x[4] >= -1);
s.t. c2: atleast 3 {i in 1..4} (x[i] >= -2);
s.t. c3: alldiff(x[3], x[2], x[4]);
s.t. c4: not ((x[4] == -2) and (x[3] != -2));
s.t. c5: abs(x[3] - x[1]) != 0;
s.t. c6: (x[1] != -1) <==> (x[4] + x[4] >= -1);
s.t. c7: not ((x[1] >= -1) and (x[2] != 0));
