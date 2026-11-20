var x{1..4} integer >= -2 <= 0;

maximize o: -4*x[1] + 3*x[2] - 3*x[3] - x[4];

s.t. c1: (x[1] < -2) and (x[1] >= -2) ==> (x[1] >= -2);
s.t. c2: not ((x[3] == -2) and (x[1] <= -2));
s.t. c3: max(x[2], x[1]) - min(x[2], x[2]) != -2;
