BUILTINS =
NPROC = 1
SH = /home/dmg/bin/sh

/home/dmg/t/%.tgz: %/changes %/xsum0.out
	tar cf - `sed "s@^\([^	]*\)	.*@$stem/\1@" $stem/xsum0.out` $stem/changes $stem/xsum0.out | gzip -c >$target
