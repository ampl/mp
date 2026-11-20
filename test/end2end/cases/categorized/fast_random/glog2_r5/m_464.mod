var x{1..5} integer >= -3 <= 2;

minimize o: 3*x[1] - 2*x[2] - 5*x[3] + 5*x[4] - 5*x[5];

s.t. c1: (x[2] >= -1) and (x[1] != 2) ==> (x[2] >= -1);
s.t. c2: abs(x[4] - x[4]) >= 0;
s.t. c3: not ((x[2] == -3) and (x[1] < 2));
s.t. c4: alldiff(x[1], x[2], x[4]);
