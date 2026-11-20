var x{1..3} integer >= 0 <= 3;

maximize o: 3*x[1] + 2*x[2] - 5*x[3];

s.t. c1: not (max(x[2], x[2]) != 1);
s.t. c2: not (max(x[2], x[2]) != 1);
s.t. c3: (if x[2] < 2 then x[2] else x[3]) < 2;
