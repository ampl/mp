var x{1..4} integer >= -2 <= 2;

minimize o: 0*x[1] + 0*x[2] - 2*x[3] - 3*x[4];

s.t. c1: ((x[4] < -2) ==> (x[4] == -1)) and ((x[4] == -1) ==> (x[1] >= -2));
s.t. c2: abs(x[3] - x[2]) == 0;
s.t. c3: atmost 2 {i in 1..4} (x[i] >= 1);
