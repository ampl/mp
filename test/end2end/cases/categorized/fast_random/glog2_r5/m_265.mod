var x{1..4} integer >= 0 <= 3;

maximize o: x[1] - 2*x[2] - x[3] - x[4];

s.t. c1: (x[1] < 0) or (x[2] <= 3) or (x[3] < 2);
s.t. c2: max(x[3], x[1]) - min(x[1], x[4]) > 3;
s.t. c3: max(x[1], x[1]) - min(x[3], x[3]) >= 0;
s.t. c4: max(x[4], x[2]) - min(x[2], x[3]) >= 2;
s.t. c5: (x[3] >= 1) <==> (x[1] + x[3] != 1);
s.t. c6: alldiff(x[1], x[2], x[3]);
