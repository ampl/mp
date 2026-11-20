var x{1..3} integer >= 0 <= 4;

minimize o: -3*x[1] - x[2] + 5*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: alldiff(x[1], x[3], x[2]);
s.t. c3: (if x[2] >= 2 then x[3] else x[2]) != 3;
s.t. c4: atmost 3 {i in 1..3} (x[i] >= 3);
s.t. c5: x[3] == 1;
s.t. c6: x[1] == 3;
