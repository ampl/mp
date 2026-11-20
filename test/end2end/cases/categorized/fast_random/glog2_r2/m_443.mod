var x{1..6} integer >= -2 <= 1;

maximize o: 5*x[1] - 4*x[2] - 2*x[3] + x[4] - 2*x[5] + 6*x[6];

s.t. c1: (if x[6] >= -2 then x[4] else x[6]) <= 0;
s.t. c2: (x[6] == 0) ==> (x[2] < -1);
s.t. c3: alldiff(x[1], x[6]);
s.t. c4: x[3] == -1;
s.t. c5: x[6] == -1;
