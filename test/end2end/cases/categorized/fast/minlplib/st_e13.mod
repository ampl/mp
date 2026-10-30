# MINLPLib instance st_e13, https://www.minlplib.org/st_e13.html
# (c) MINLPLib, licensed under CC BY 4.0,
# https://creativecommons.org/licenses/by/4.0/ ; only this header added.
#
#  MINLP written by GAMS Convert at 01/12/18 13:43:22
#  
#  Equation counts
#      Total        E        G        L        N        X        C        B
#          3        1        0        2        0        0        0        0
#  
#  Variable counts
#                   x        b        i      s1s      s2s       sc       si
#      Total     cont   binary  integer     sos1     sos2    scont     sint
#          3        2        1        0        0        0        0        0
#  FX      0        0        0        0        0        0        0        0
#  
#  Nonzero counts
#      Total    const       NL      DLL
#          7        6        1        0
# 
#  Reformulation has removed 1 variable and 1 equation


var b1 binary >= 0, <= 1;
var x2 >= 0, <= 1.6;

minimize obj:    b1 + 2*x2;

subject to

e1: -x2^2 - b1 <= -1.25;

e2:    b1 + x2 <= 1.6;
