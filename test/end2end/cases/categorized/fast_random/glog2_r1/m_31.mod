var x{1..4} integer >= -2 <= 2;

minimize o: 4*x[1] - 5*x[2] - x[3] + 3*x[4];

s.t. c1: ((x[1] > -2) ==> (x[1] != 2)) and ((x[1] != 2) ==> (x[4] >= -2));
s.t. c2: not ((x[4] < -2) and (x[3] <= -3));
s.t. c3: exists {i in 1..4} x[i] == 2;
s.t. c4: (x[3] <= -3) ==> (x[1] > 2) else (x[3] < 2);
