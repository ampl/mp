# MINLPLib instance synthes1, https://www.minlplib.org/synthes1.html
# (c) MINLPLib, licensed under CC BY 4.0,
# https://creativecommons.org/licenses/by/4.0/ ; only this header added.
#
#  MINLP written by GAMS Convert at 01/12/18 13:43:53
#  
#  Equation counts
#      Total        E        G        L        N        X        C        B
#          7        1        2        4        0        0        0        0
#  
#  Variable counts
#                   x        b        i      s1s      s2s       sc       si
#      Total     cont   binary  integer     sos1     sos2    scont     sint
#          7        4        3        0        0        0        0        0
#  FX      0        0        0        0        0        0        0        0
#  
#  Nonzero counts
#      Total    const       NL      DLL
#         23       17        6        0
# 
#  Reformulation has removed 1 variable and 1 equation


var x1 >= 0, <= 2;
var x2 >= 0, <= 2;
var x3 >= 0, <= 1;
var b4 binary >= 0, <= 1;
var b5 binary >= 0, <= 1;
var b6 binary >= 0, <= 1;

minimize obj: (-18*log(1 + x2)) - 19.2*log(1 + x1 - x2) + 10*x1 - 7*x3 + 5*b4
     + 6*b5 + 8*b6 + 10;

subject to

e2: 0.8*log(1 + x2) + 0.96*log(1 + x1 - x2) - 0.8*x3 >= 0;

e3: log(1 + x2) + 1.2*log(1 + x1 - x2) - x3 - 2*b6 >= -2;

e4:  - x1 + x2 <= 0;

e5:    x2 - 2*b4 <= 0;

e6:    x1 - x2 - 2*b5 <= 0;

e7:    b4 + b5 <= 1;
