var x{1..5} integer >= 0 <= 3;

maximize o: -4*x[1] - 5*x[2] + 3*x[3] + 4*x[4] + 5*x[5];

s.t. c1: (x[4] >= 3) and ((x[1] > 3) or not (x[4] >= 3));
s.t. c2: not (x[3] == 3);
s.t. c3: (x[3] == 3) ==> ((x[4] >= 3) or (x[1] > 3));
