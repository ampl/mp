var x{1..3} integer >= 0 <= 2;

maximize o: 2*x[1] + 2*x[2] + 6*x[3];

s.t. c1: exists {i in 1..3} x[i] == 1;
s.t. c2: alldiff(x[2], x[1]);
s.t. c3: (x[3] != 2) or (x[3] < 0) or (x[3] != 2);
s.t. c4: (if x[1] < 1 then x[2] else x[1]) < 2;
s.t. c5: x[2] == 1;
