var x{1..3} integer >= 0 <= 4;

maximize o: -5*x[1] - 5*x[2] + 5*x[3];

s.t. c1: (x[3] <= 4) <==> (x[2] + x[2] < 4);
s.t. c2: numberof 0 in (x[1], x[2], x[3]) != 1;
s.t. c3: alldiff(x[2], x[1]);
