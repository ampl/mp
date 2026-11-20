var x{1..3} integer >= 0 <= 3;

maximize o: 2*x[1] + 3*x[2] - 2*x[3];

s.t. c1: (x[2] < 3) and (x[2] != 2) ==> (x[3] >= 3);
s.t. c2: abs(x[1] - x[3]) > 3;
s.t. c3: not ((x[1] >= 3) and (x[3] != 0));
s.t. c4: x[1] == 3;
