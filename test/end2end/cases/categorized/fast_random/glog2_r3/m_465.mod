var x{1..5} integer >= 0 <= 3;

minimize o: -4*x[1] - 2*x[2] - 5*x[3] - x[4] + x[5];

s.t. c1: not ((x[4] <= 3) and (x[1] <= 3));
s.t. c2: atmost 2 {i in 1..5} (x[i] >= 2);
s.t. c3: not ((x[2] >= 3) and (x[2] <= 3));
s.t. c4: abs(x[3] - x[2]) <= 0;
s.t. c5: (if x[2] > 1 then x[2] else x[2]) > 2;
s.t. c6: (if x[1] != 2 then x[4] else x[5]) == 2;
s.t. c7: atleast 3 {i in 1..5} (x[i] >= 0);
