var x{1..5} integer >= 0 <= 4;

maximize o: -2*x[1] - x[2] - 4*x[3] + 4*x[4] - 2*x[5];

s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) != 0;
s.t. c2: atleast 1 {i in 1..5} (x[i] >= 3);
s.t. c3: max(x[5], x[3]) - min(x[5], x[5]) <= 2;
s.t. c4: ((x[1] >= 0) ==> (x[1] == 3)) and ((x[1] == 3) ==> (x[2] >= 0));
s.t. c5: alldiff(x[3], x[4]);
s.t. c6: alldiff(x[4], x[2]);
