var x{1..5} integer >= 0 <= 3;

minimize o: 6*x[1] - x[2] - x[3] - 3*x[4] + x[5];

s.t. c1: (x[1] != 1) ==> (x[5] != 3) else (x[1] >= 0);
s.t. c2: ((x[3] == 2) ==> (x[3] < 0)) and ((x[3] < 0) ==> (x[3] >= 2));
s.t. c3: exists {i in 1..5} x[i] == 1;
s.t. c4: (x[1] != 1) ==> (x[5] != 3) else (x[1] >= 0);
s.t. c5: exists {i in 1..5} x[i] == 1;
