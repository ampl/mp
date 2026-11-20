var x{1..5} integer >= -2 <= 3;

maximize o: 5*x[1] + 5*x[2] - 2*x[3] + 0*x[4] - 5*x[5];

s.t. c1: (if x[4] > -2 then x[4] else x[2]) == 3;
s.t. c2: max(x[4], x[1]) - min(x[3], x[4]) == -1;
s.t. c3: alldiff(x[4], x[3], x[5], x[1]);
s.t. c4: ((x[2] <= 3) ==> (x[2] != 0)) and ((x[2] != 0) ==> (x[3] >= 3));
