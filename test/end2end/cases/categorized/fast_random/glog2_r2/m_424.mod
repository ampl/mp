var x{1..5} integer >= -3 <= 2;

maximize o: -2*x[1] + 3*x[2] + 2*x[3] + 3*x[4] + 3*x[5];

s.t. c1: (x[3] == -3) <==> (x[1] + x[2] <= 0);
s.t. c2: alldiff(x[3], x[2], x[1], x[5]);
s.t. c3: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) < 1;
s.t. c4: (x[3] < -2) ==> (x[2] >= 2);
