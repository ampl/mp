var x{1..3} integer >= 0 <= 4;

minimize o: 6*x[1] - x[2] + 2*x[3];

s.t. c1: (x[2] > 2) and (x[1] >= 0) ==> (x[1] >= 2);
s.t. c2: (x[1] <= 3) or (x[1] == 1) or (x[3] != 1);
s.t. c3: (if x[2] >= 0 then x[2] else x[2]) >= 0;
s.t. c4: (x[3] > 3) and (x[3] < 2) ==> (x[3] >= 3);
s.t. c5: not ((x[2] == 2) and (x[1] == 0));
