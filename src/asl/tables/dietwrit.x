# Write, then read files {calories, foods, nutrient}.{tab, bit, abt}
# Corresponding .bit and .abt files should have the same contents.

model diet.mod;

set Ttype;
set Handlers;

table calories{t in Ttype} ('calories.' & t): [NUTR,FOOD] amt;

table foods{t in Ttype} ('foods.' & t):
	[FOOD] INOUT cost, f_min, f_max;

table nutrients{t in Ttype} ('nutrients.' & t):
	[NUTR] INOUT, n_min, n_max;

data diet2a.dat;
data; set Ttype := tab bit abt;
set Handlers := simpbit.dll fullbit.dll;

load simpbit.dll;
write table {t in Ttype} calories[t];
write table {t in Ttype} foods[t];
write table {t in Ttype} nutrients[t];
unload simpbit.dll;

for {h in Handlers} {
	load (h);
	for{t in Ttype} {
		if h == 'fullbit.dll' && t == 'abt' then break;
		display h, t;
		reset data amt, cost, f_min, f_max, n_min, n_max;
		read table calories[t];	# will abort on (fullbit.dll, abt)
		read table foods[t];
		read table nutrients[t];
		solve;
		display Buy;
		}
	unload (h);
	}
