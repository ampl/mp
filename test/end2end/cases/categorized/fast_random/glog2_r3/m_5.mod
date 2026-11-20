var x{1..3} integer >= 0 <= 2;

maximize o: -2*x[1] - 3*x[2] + 4*x[3];

s.t. c1: max(x[1], x[2]) - min(x[3], x[2]) != 0;
s.t. c2: (x[3] >= 0) ==> (x[1] > 2);
s.t. c3: max(x[2], x[3]) - min(x[2], x[1]) != 0;
s.t. c4: abs(x[2] - x[2]) == 1;
s.t. c5: (if x[2] <= 2 then x[2] else x[1]) != 1;
s.t. c6: (x[3] < 0) ==> (x[3] < 1);
s.t. c7: abs(x[2] - x[3]) > 0;
