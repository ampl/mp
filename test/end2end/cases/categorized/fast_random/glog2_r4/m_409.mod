var x{1..5} integer >= 0 <= 2;
minimize o: (1)*x[1] + (3)*x[2] + (4)*x[3] + (1)*x[4] + (-3)*x[5];
s.t. c1: not ((x[3] <= 0) and (x[2] == 1));
s.t. c2: ((x[2] < 0) ==> (x[4] >= 0)) and ((x[4] >= 0) ==> (x[5] >= 0));
s.t. c3: alldiff(x[5], x[1]);
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 2;
s.t. c5: not ((x[3] <= 0) and (x[2] == 1));
s.t. c6: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 2;
