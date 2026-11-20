var x{1..5} integer >= 0 <= 4;

maximize o: -x[1] + 5*x[2] - 2*x[3] - 5*x[4] + 0*x[5];

s.t. c1: (x[4] == 3) <==> (x[5] + x[1] >= 4);
s.t. c2: ((x[3] > 2) ==> (x[5] != 1)) and ((x[5] != 1) ==> (x[1] >= 2));
s.t. c3: not ((x[5] == 4) and (x[5] >= 3));
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
s.t. c5: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) >= 3;
s.t. c6: alldiff(x[4], x[3], x[2]);
s.t. c7: (if x[2] < 4 then x[3] else x[5]) < 3;
