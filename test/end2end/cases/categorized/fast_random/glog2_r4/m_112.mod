var x{1..5} integer >= 0 <= 3;

maximize o: -3*x[1] + 3*x[2] + x[3] - 5*x[4] + 5*x[5];

s.t. c1: alldiff(x[1], x[4], x[5]);
s.t. c2: max(x[4], x[4]) - min(x[4], x[5]) != 0;
s.t. c3: max(x[4], x[1]) - min(x[1], x[1]) > 1;
s.t. c4: (x[3] != 3) ==> (x[2] == 2) else (x[2] > 2);
