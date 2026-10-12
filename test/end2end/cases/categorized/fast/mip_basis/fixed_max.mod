# Exercises a negative integer, maximization, an objective constant,
# an equality, and a ranged row. The fixed LP range dual is 3.
var k integer >= -3, <= -1;
var b binary;
var y >= 0, <= 4;
maximize profit: 5 + 4*k + b + 3*y;
subject to range: 1 <= y + k <= 2;
subject to binary_value: b = 1;
option presolve 0;
