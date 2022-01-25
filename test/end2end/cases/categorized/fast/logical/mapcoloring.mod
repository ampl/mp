param NumColors;

set Countries;
set Neighbors within Countries cross Countries;

var color{Countries} integer >= 1 <= NumColors;

s.t. different_colors{(c1, c2) in Neighbors}:
  color[c1] != color[c2];

data;

param NumColors := 4;

set Countries := Belgium Denmark France Germany Luxembourg Netherlands;

set Neighbors :=
  Belgium France 
  Belgium Germany 
  Belgium Netherlands
  Belgium Luxembourg
  Denmark Germany 
  France  Germany 
  France  Luxembourg
  Germany Luxembourg
  Germany Netherlands;
