# MINLPLib instance ex1223, https://www.minlplib.org/ex1223.html
# (c) MINLPLib, licensed under CC BY 4.0,
# https://creativecommons.org/licenses/by/4.0/ ; only this header added.
#
#  MINLP written by GAMS Convert at 01/12/18 13:31:49
#  
#  Equation counts
#      Total        E        G        L        N        X        C        B
#         14        5        0        9        0        0        0        0
#  
#  Variable counts
#                   x        b        i      s1s      s2s       sc       si
#      Total     cont   binary  integer     sos1     sos2    scont     sint
#         12        8        4        0        0        0        0        0
#  FX      0        0        0        0        0        0        0        0
#  
#  Nonzero counts
#      Total    const       NL      DLL
#         40       23       17        0
# 
#  Reformulation has removed 1 variable and 1 equation


var x1 >= 0, <= 10;
var x2 >= 0, <= 10;
var x3 >= 0, <= 10;
var x4 >= 0, <= 1;
var x5 >= 0, <= 1;
var x6 >= 0, <= 1;
var x7 >= 0, <= 1;
var b8 binary >= 0, <= 1;
var b9 binary >= 0, <= 1;
var b10 binary >= 0, <= 1;
var b11 binary >= 0, <= 1;

minimize obj: (-1 + x4)^2 + (-2 + x5)^2 + (-1 + x6)^2 - log(1 + x7) + (-1 + x1)
    ^2 + (-2 + x2)^2 + (-3 + x3)^2;

subject to

e1:    x1 + x2 + x3 + b8 + b9 + b10 <= 5;

e2: x6^2 + x1^2 + x2^2 + x3^2 <= 5.5;

e3:    x1 + b8 <= 1.2;

e4:    x2 + b9 <= 1.8;

e5:    x3 + b10 <= 2.5;

e6:    x1 + b11 <= 1.2;

e7: x5^2 + x2^2 <= 1.64;

e8: x6^2 + x3^2 <= 4.25;

e9: x5^2 + x3^2 <= 4.64;

e10:    x4 - b8 = 0;

e11:    x5 - b9 = 0;

e12:    x6 - b10 = 0;

e13:    x7 - b11 = 0;
