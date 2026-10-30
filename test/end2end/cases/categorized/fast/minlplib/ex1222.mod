# MINLPLib instance ex1222, https://www.minlplib.org/ex1222.html
# (c) MINLPLib, licensed under CC BY 4.0,
# https://creativecommons.org/licenses/by/4.0/ ; only this header added.
#
#  MINLP written by GAMS Convert at 01/12/18 13:31:48
#  
#  Equation counts
#      Total        E        G        L        N        X        C        B
#          4        1        0        3        0        0        0        0
#  
#  Variable counts
#                   x        b        i      s1s      s2s       sc       si
#      Total     cont   binary  integer     sos1     sos2    scont     sint
#          4        3        1        0        0        0        0        0
#  FX      0        0        0        0        0        0        0        0
#  
#  Nonzero counts
#      Total    const       NL      DLL
#          9        7        2        0
# 
#  Reformulation has removed 1 variable and 1 equation


var x1 := 0.2, >= 0.2, <= 1;
var x2 := -1, >= -2.22554, <= -1;
var b3 binary >= 0, <= 1;

minimize obj: 5*(-0.5 + x1)^2 - 0.7*b3 + 0.8;

subject to

e1: -exp((-0.2) + x1) - x2 <= 0;

e2:    x2 + 1.1*b3 <= -1;

e3:    x1 - 1.2*b3 <= 0;
