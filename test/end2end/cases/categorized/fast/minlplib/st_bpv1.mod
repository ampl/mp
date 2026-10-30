# MINLPLib instance st_bpv1, https://www.minlplib.org/st_bpv1.html
# (c) MINLPLib, licensed under CC BY 4.0,
# https://creativecommons.org/licenses/by/4.0/ ; only this header added.
#
#  NLP written by GAMS Convert at 01/12/18 13:43:19
#  
#  Equation counts
#      Total        E        G        L        N        X        C        B
#          5        1        3        1        0        0        0        0
#  
#  Variable counts
#                   x        b        i      s1s      s2s       sc       si
#      Total     cont   binary  integer     sos1     sos2    scont     sint
#          5        5        0        0        0        0        0        0
#  FX      0        0        0        0        0        0        0        0
#  
#  Nonzero counts
#      Total    const       NL      DLL
#         13        9        4        0
# 
#  Reformulation has removed 1 variable and 1 equation


var x1 >= 0, <= 27;
var x2 >= 0, <= 16;
var x3 >= 0, <= 10;
var x4 >= 0, <= 10;

minimize obj: x1*x3 + x2*x4;

subject to

e1:    x1 + 3*x2 >= 30;

e2:    2*x1 + x2 >= 20;

e3:  - 1.6667*x3 + x4 >= 10;

e4:    x3 + x4 <= 15;
