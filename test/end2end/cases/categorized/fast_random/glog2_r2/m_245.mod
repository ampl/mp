var x{1..3} integer >= -2 <= 1;

minimize o: -3*x[1] + x[2] - 4*x[3];

s.t. c1: (x[2] <= 0) and (x[2] >= 1) ==> (x[2] >= 0);
s.t. c2: (x[2] >= 1) ==> (x[1] == 1);
s.t. c3: (if x[3] >= -2 then x[1] else x[3]) > -2;
s.t. c4: exists {i in 1..3} x[i] == -1;
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: atleast 1 {i in 1..3} (x[i] >= -2);
s.t. c7: ((x[2] >= 0) ==> (x[1] < -2)) and ((x[1] < -2) ==> (x[1] >= 0));
