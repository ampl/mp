var x{1..5} integer >= 0 <= 4;

minimize o: 2*x[1] + 2*x[2] + 5*x[3] - 4*x[4] - x[5];

s.t. c1: (x[4] == 1) <==> (x[1] + x[4] == 2);
s.t. c2: ((x[5] > 0) ==> (x[1] != 4)) and ((x[1] != 4) ==> (x[3] >= 0));
s.t. c3: exists {i in 1..5} x[i] == 0;
s.t. c4: (if x[2] == 2 then x[3] else x[4]) != 2;
s.t. c5: (if x[4] < 0 then x[1] else x[4]) <= 1;
s.t. c6: not ((x[3] == 3) and (x[5] > 1));
s.t. c7: forall {i in 1..4} x[i] != x[i+1];
