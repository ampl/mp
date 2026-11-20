var x{1..3} integer >= 0 <= 3;

minimize o: -5*x[1] + 2*x[2] + x[3];

s.t. c1: ((x[2] == 3) ==> (x[3] >= 0)) and ((x[3] >= 0) ==> (x[1] >= 3));
s.t. c2: numberof 3 in (x[1], x[2], x[3]) == 3;
s.t. c3: atleast 2 {i in 1..3} (x[i] >= 4);
s.t. c4: (x[3] > 3) and (x[3] <= 4) ==> (x[3] >= 3);
