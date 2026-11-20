# The quadratic objective term x^2 has a coefficient other than 1 and
# shares x^2 with a quadratic constraint (recognized as a cone): the
# objective's rotated cone must use sqrt(|coefficient|), not the coefficient.
# Optimum: x = 3/(2*0.8) = 1.875, y = 0, objective 2.8125.

var x >= 0 <= 10;
var y >= 0 <= 10;
maximize qobj: 3*x - 0.8*x^2 - y;
subj to qconstr: x^2 + y^2 <= 100;
