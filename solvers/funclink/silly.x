# Demonstration of functions linked into AMPL.

function ginv;		# ginv(x) returns the generalized inverse of x:
			# ginv(x) = if x == 0 then 0 else 1/x

display{i in -3..3} ginv(i);


function sginv symbolic; # symbolic (character-valued) variant of ginv

display{i in -3..3} sginv(i);


function hypot(Reals,Reals); # hypot(x,y) = sqrt(x^2 + y^2);

display hypot(3,4);


function ncall();		# returns its invocation count

display {1..3} ncall();


function rncall() random;	# returns its invocation count

display {1..3} rncall();

# ncall() and rncall() are linked to the same function;
# but AMPL thinks ncall() always returns the same value,
# because it is not declared random.  Thus display{1..3} ncall()
# prints 1 three times, and a subsequent display{1..3} rncall()
# prints 2 3 4.


function mean0;		# mean of its args, which must be numeric

display mean0(2,3,4);


function mean;		# mean of its args, with symbolic args
			# treated as 0 (after a complaint)

display mean('a',2,3,4);


function kth symbolic;	# kth(k,a1,a2,...,an) returns ak

display{i in 1..3} kth(i,'a',i+10,'Last arg');

function ginvae;	# like ginv, but enrolls functions
			# that are called at "reset;" and the
			# end of execution

display{i in 0..2} ginvae(i);

reset;			# causes "Got to At_reset" messages

function ginvae;

display {i in -1 .. 1} ginvae(i);

quit;			# causes "Got to At_reset" and "Got to At_end"
			# messages
