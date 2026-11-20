var x{1..5} integer >= 0 <= 3;

maximize o: 3000*x[1] + 51*x[2] + 3000*x[3] - 666*x[4] + 1000*x[5];

s.t. c1: abs(x[3] - x[3]) >= 3;
s.t. c2: abs(x[3] - x[1]) <= 1;
