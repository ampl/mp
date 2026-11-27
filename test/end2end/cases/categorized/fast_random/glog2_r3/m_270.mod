var x{1..5} integer >= 0 <= 4;

minimize o: x[1] - 4*x[2] + 3*x[3] + 5*x[4] + 4*x[5];

s.t. c1: alldiff(x[5], x[2], x[3]);
s.t. c2: (if x[4] != 0 then x[2] else x[5]) <= 0;
s.t. c3: ((x[4] <= 0) ==> (x[4] < 4)) and ((x[4] < 4) ==> (x[5] >= 0));
s.t. c4: (if x[5] == 2 then x[1] else x[1]) < 3;
