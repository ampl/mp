var x{1..4} integer >= 0 <= 4;

minimize o: 0*x[1] + 5*x[2] + 0*x[3] + 0*x[4];

s.t. c1: ((x[4] > 3) ==> (x[2] == 4)) and ((x[2] == 4) ==> (x[3] >= 3));
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: abs(x[4] - x[1]) < 3;
s.t. c4: (if x[2] < 0 then x[3] else x[3]) < 4;
s.t. c5: ((x[4] > 3) ==> (x[2] == 4)) and ((x[2] == 4) ==> (x[3] >= 3));
s.t. c6: (if x[2] < 0 then x[3] else x[3]) < 4;
