var x{1..4} integer >= 0 <= 4;

maximize o: 2*x[1] - 4*x[2] - x[3] + 4*x[4];

s.t. c1: (x[1] < 4) and (x[4] < 3) ==> (x[1] >= 4);
s.t. c2: not ((x[3] > 3) and (x[2] < 4));
s.t. c3: (if x[4] != 4 then x[1] else x[1]) < 3;
s.t. c4: (x[1] >= 4) ==> (x[2] < 2);
s.t. c5: forall {i in 1..3} x[i] != x[i+1];
s.t. c6: (x[3] != 4) ==> (x[3] == 4) else (x[4] > 4);
