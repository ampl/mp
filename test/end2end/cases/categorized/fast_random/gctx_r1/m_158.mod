var x{1..3} integer >= -1 <= 1;

minimize o: 3*x[1] + 5*x[2] + 4*x[3];

s.t. c1: (x[2] + x[1] != 0) ==> (x[3] + x[3] < 1);
s.t. c2: not ((x[2] < x[2]) and (x[2] + x[1] != 0));
s.t. c3: (if x[2] + x[1] != 0 then x[1] else x[3]) < 0;
