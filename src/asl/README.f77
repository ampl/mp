The following changes (to solvers/makefile and solvers/*/makefile)
may permit use of the native Fortran 77 compiler on some systems.
On other systems, you may be able to discover suitable changes by
studying system documents (including man pages), and by using the nm
command to look at the names in relevant libraries.  Often the
compiler flags "-v -v" (two -v's) will cause the compiler to tell
you what libraries it references; if, say, /usr/lib/libF77.a is
among them, you could try thing like

	nm /usr/lib/libF77.a | grep -i arg

to get a hint at the system's variant (if any) of xargv.

If you make the changes shown below, also remove -lf2c from
solvers/*/makefile.  In other words, if you use the native Fortran
compiler, do not link against libf2c.a.


Sun SunOS:
	CFLAGS = -O -DKR_headers -DMAIN__=MAIN_ -Dxargv=_xargv


Sun Solaris:
	CFLAGS = -O -DMAIN__=main_ -Dxargv=__xargv
	For solvers defining MAIN__, add fmain.o built
	from fmain.f consisting of the two lines
		call main
		end


HP:
	CFLAGS = -Aa -O
	FFLAGS = +ppu
	For solvers defining MAIN__, add fmain.o built
	from fmain.c consisting of the following:

		char **xargv;
		extern void MAIN__(void);
		main(int argc, char **argv)
		{
			xargv = argv;
			MAIN__();
			return 0;
			}


IBM RS6000:
	CFLAGS = -O -Dxargv=p_xargv -DMAIN__=main_
	FFLAGS = -qextname
	For solvers defining MAIN__, add fmain.o built
	from fmain.f consisting of the two lines
		call main
		end


SGI IRIX:
	CFLAGS = -O -Dxargv=f77argv


DEC OSF1 (Unix for Alpha chip):
	CFLAGS = -O -Dxargv=for__a_argv

Linux (with g77):
	CFLAGS = -O -Dxargv=f__xargv

Linux (with gfortran):
	CFLAGS = -O
	For solvers defining MAIN__, supply fmain.o as for HP above.
	For solvers (e.g., snopt) that reference etime_(), append

		extern double xectim_();
		float etime_(float *tarray) { return (float) xectim_(); }

	to fmain.c (source for fmain.o; see HP above).
