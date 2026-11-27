var x{1..3} integer >= -2 <= 1;

maximize o: -333*x[1] + 2000*x[2] + 3000*x[3];

s.t. c1: abs(x[1] - x[2]) == 3;
s.t. c2: (x[3] < -1) and (x[3] >= 1) ==> (x[3] >= -1);
