var x{1..6} integer >= 0 <= 3;

maximize o: 2*x[1] - x[2] - 5*x[3] - x[4] + 6*x[5] - x[6];

s.t. c1: (x[4] >= 2) ==> (x[3] == 3) else (x[1] != 2);
s.t. c2: (if x[5] == 1 then x[4] else x[2]) < 0;
s.t. c3: alldiff(x[3], x[6]);
s.t. c4: ((x[6] >= 0) ==> (x[2] == 0)) and ((x[2] == 0) ==> (x[3] >= 0));
