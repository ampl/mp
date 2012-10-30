# A constraint programming model to solve the SEND + MORE = MONEY puzzle.

set Letters;
var d{Letters} >= 0 <= 9 integer;

s.t. nonzeroS: d['S'] != 0;
s.t. nonzeroM: d['M'] != 0;

s.t. equation:     1000 * d['S'] + 100 * d['E'] + 10 * d['N'] + d['D'] +
                   1000 * d['M'] + 100 * d['O'] + 10 * d['R'] + d['E'] =
  10000 * d['M'] + 1000 * d['O'] + 100 * d['N'] + 10 * d['E'] + d['Y'];

s.t. different: alldiff{l in Letters} d[l];

data;

set Letters := S E N D M O R Y;
