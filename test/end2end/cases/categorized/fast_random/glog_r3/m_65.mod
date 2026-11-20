var x{1..5} integer >= 0 <= 4;

maximize o: 5*x[1] - 3*x[2] + 2*x[3] + 5*x[4] + x[5];

s.t. c1: (x[2] >= 2) ==> (x[2] != 4);
s.t. c2: ((x[1] != 2) ==> (x[5] == 2)) and ((x[5] == 2) ==> (x[5] >= 2));
s.t. c3: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) == 3;
s.t. c4: alldiff(x[1], x[2], x[5]);
