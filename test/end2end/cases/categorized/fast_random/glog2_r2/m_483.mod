var x{1..3} integer >= 0 <= 4;

minimize o: 6*x[1] + 6*x[2] + 5*x[3];

s.t. c1: (if x[1] == 2 then x[1] else x[2]) <= 4;
s.t. c2: (x[3] == 1) and (x[2] < 1) ==> (x[2] >= 1);
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: ((x[1] < 3) ==> (x[1] == 2)) and ((x[1] == 2) ==> (x[3] >= 3));
s.t. c5: x[1] == 3;
