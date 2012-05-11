## Makefile for WATCOM C/C++ 11.0a (for making Windows 9x and NT .exe files)
## amplfunc.dll provides user-defined functions.

# Invoke with "wmake -u -f makefile.wat" .

# $S = ampl/solvers directory
S = ..
CC = wcc386
CFLAGS = -I$S -fpd

.c.obj:
	$(CC) $(CFLAGS) $*.c

## Sample amplfunc.dll (which also creates amplfunc.lib):
amplfunc.dll: libmain.obj funcadd.obj amplfunc.lnk
	wlink @amplfunc
	ren funcadd.dll amplfunc.dll

## Sample solver creation...

# $(myobjects) = list of .obj files
#myobjects = ....

mysolver.exe: $(myobjects)
	wcl386 -l=nt -bt=nt -fe=mysolver.exe $(myobjects) $S\amplsolv.lib
