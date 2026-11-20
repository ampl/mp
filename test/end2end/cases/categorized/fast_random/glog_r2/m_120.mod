var x{1..5} integer >= -2 <= 2;

minimize o: 0*x[1] - 2*x[2] - 5*x[3] - 4*x[4] + 0*x[5];

s.t. c1: abs(x[2] - x[5]) != 3;
s.t. c2: (x[3] <= 0) ==> (x[2] < -2) else (x[3] == 1);
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) < 2;
s.t. c4: alldiff(x[1], x[4], x[2]);
