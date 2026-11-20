var x{1..6} integer >= 0 <= 3;

maximize o: 0*x[1] - x[2] + 5*x[3] + 5*x[4] + 5*x[5] + 3*x[6];

s.t. c1: alldiff(x[4], x[1], x[3], x[2]);
s.t. c2: max(x[6], x[1]) - min(x[5], x[6]) < 3;
s.t. c3: (x[2] == 0) ==> (x[5] <= 0);
s.t. c4: x[5] == 1;
s.t. c5: x[1] == 3;
