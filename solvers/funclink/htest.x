var x{i in 1..2} := i;
function hypot;
minimize zot: (hypot({i in 1..2} x[i]) - 5)^2 + (x[1] - .2)^2;
solve;
display x;
print zot;
