var x{1..5} integer >= 0 <= 4;

minimize o: -2000*x[1] + 2000*x[2] + 2000*x[3] - 666*x[4] + 3000*x[5];

s.t. c1: (if x[2] <= 3 then x[1] else x[3]) >= 3;
s.t. c2: alldiff(x[4], x[5], x[1], x[2]);
s.t. c3: (if x[1] <= 1 then x[3] else x[2]) != 2;
