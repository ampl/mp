var x{1..5} integer >= -2 <= 0;

maximize o: 5*x[1] + 5*x[2] + 3*x[3] + 0*x[4] + 3*x[5];

s.t. c1: (x[4] <= -3) ==> (x[3] > 0);
s.t. c2: max(x[5], x[1]) - min(x[2], x[1]) != 0;
s.t. c3: numberof -3 in (x[1], x[2], x[3], x[4], x[5]) >= 0;
