var x{1..3} integer >= 0 <= 4;

maximize o: 4*x[1] - x[2] + 5*x[3];

s.t. c1: abs(x[3] - x[1]) > 1;
s.t. c2: (x[3] <= 1) <==> (x[2] + x[3] >= 1);
s.t. c3: alldiff(x[3], x[1]);
s.t. c4: (if x[1] < 0 then x[2] else x[3]) <= 3;
