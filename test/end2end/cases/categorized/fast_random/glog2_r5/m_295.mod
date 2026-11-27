var x{1..5} integer >= 0 <= 3;

maximize o: -4*x[1] + 3*x[2] + 0*x[3] + 0*x[4] - 5*x[5];

s.t. c1: (x[2] < 0) and (x[5] >= 1) ==> (x[4] >= 0);
s.t. c2: ((x[4] != 2) ==> (x[4] >= 1)) and ((x[4] >= 1) ==> (x[4] >= 2));
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) != 2;
s.t. c5: (x[1] < 2) or (x[3] != 3) or (x[4] == 1);
s.t. c6: alldiff(x[5], x[2]);
