var x{1..5} integer >= -2 <= 2;

maximize o: 5*x[1] - 4*x[2] - 2*x[3] + 4*x[4] + 0*x[5];

s.t. c1: (x[2] == -2) or (x[3] > 1) or (x[5] < -1);
s.t. c2: abs(x[5] - x[2]) >= 4;
