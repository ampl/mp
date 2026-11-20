var x{1..3} integer >= 0 <= 2;

maximize o: -666*x[1] - 1000*x[2] + 1000*x[3];

s.t. c1: abs(x[2] - x[1]) >= 2;
s.t. c2: max(x[1], x[3]) - min(x[1], x[2]) != 0;
s.t. c3: (x[3] != 0) and (x[3] == 1) ==> (x[1] >= 0);
s.t. c4: (if x[1] >= 1 then x[3] else x[3]) != 2;
