var x{1..3} integer >= 0 <= 4;

minimize o: 3*x[1] - 2*x[2] - 5*x[3];

s.t. c1: (x[2] < 4) or (x[1] <= 3) or (x[2] == 1);
s.t. c2: ((x[3] == 4) ==> (x[3] == 2)) and ((x[3] == 2) ==> (x[2] >= 4));
s.t. c3: exists {i in 1..3} x[i] == 1;
s.t. c4: numberof 2 in (x[1], x[2], x[3]) != 0;
s.t. c5: x[2] == 1;
s.t. c6: x[3] == 3;
