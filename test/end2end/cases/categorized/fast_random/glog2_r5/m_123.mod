var x{1..5} integer >= -2 <= 2;

maximize o: 3*x[1] - x[2] - 4*x[3] - 2*x[4] - 5*x[5];

s.t. c1: alldiff(x[2], x[3], x[1], x[4]);
s.t. c2: abs(x[4] - x[5]) < 1;
s.t. c3: alldiff(x[3], x[2], x[4]);
s.t. c4: x[5] == 2;
s.t. c5: x[3] == 0;
