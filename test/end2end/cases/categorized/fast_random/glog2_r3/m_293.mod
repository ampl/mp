var x{1..5} integer >= -2 <= 2;

maximize o: 6*x[1] - 2*x[2] + 4*x[3] - 2*x[4] + 0*x[5];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) != 0;
s.t. c2: not ((x[3] <= 2) and (x[4] < -2));
s.t. c3: (x[5] < 2) and (x[4] >= 1) ==> (x[1] >= 2);
s.t. c4: x[5] == 1;
s.t. c5: x[1] == 2;
