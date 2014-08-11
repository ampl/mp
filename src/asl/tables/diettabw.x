# write calories.tab, foods.tab, nutrients.tab.

model diet.mod;

table calories: [NUTR,FOOD] amt;

table foods:
	[FOOD] INOUT cost, f_min, f_max;

table nutrients:
	[NUTR] INOUT, n_min, n_max;

data diet2a.dat;

write table calories;
write table foods;
write table nutrients;
