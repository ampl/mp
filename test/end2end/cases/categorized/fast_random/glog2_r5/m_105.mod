var x{1..3} integer >= -2 <= 0;

minimize o: 2*x[1] + 0*x[2] + 4*x[3];

s.t. c1: numberof 0 in (x[1], x[2], x[3]) == 3;
s.t. c2: exists {i in 1..3} x[i] == -2;
s.t. c3: ((x[2] < 0) ==> (x[2] <= 0)) and ((x[2] <= 0) ==> (x[2] >= 0));
s.t. c4: atleast 1 {i in 1..3} (x[i] >= -2);
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: ((x[2] <= -1) ==> (x[2] >= -1)) and ((x[2] >= -1) ==> (x[1] >= -1));
