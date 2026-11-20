var x{1..7} integer >= 0 <= 1;

minimize o: 5*x[1] + 5*x[2] + 0*x[3] + 3*x[4] + 5*x[5] - 5*x[6] - 4*x[7];

s.t. c1: atmost 2 {i in 1..7} (x[i] >= 1);
s.t. c2: (x[7] >= 0) or (x[4] < 1) or (x[4] > 0);
s.t. c3: alldiff(x[4], x[6]);
s.t. c4: ((x[4] >= 1) ==> (x[2] == 0)) and ((x[2] == 0) ==> (x[5] >= 1));
