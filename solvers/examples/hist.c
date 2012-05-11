/****************************************************************
Copyright (C) 2002 AMPL Optimization LLC
All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all source copies
and that both that the copyright notice and this permission notice and
warranty disclaimer appear in supporting documentation, and that the
name of AMPL Optimization LLC or any of its entities not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

AMPL Optimization LLC DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS.  IN NO EVENT SHALL AMPL Optimization LLC OR ANY OF ITS
ENTITIES BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES
OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
SOFTWARE.

Author: David M. Gay (dmg at ampl dot com)

****************************************************************/

/* Supply bash-style history to any program that reads stdin. */

/* This uses the GNU history and readline libraries: link with */
/* -lhistory -lreadline -lcurses */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <termios.h>
#include <sys/ioctl.h>

 static char *progname;
 static int in_readline, pipegone;

 static void
catch_sigpipe(int sig)
{
	if (!pipegone) {
		pipegone = 1;
		if (in_readline)
			ioctl(0, TIOCSTI, "\n"); /* force readline() to complete */
		}
	}

 static int
usage(int rc)
{
	FILE *f = rc ? stderr : stdout;
	fprintf(f, "Usage: %s [-nnn] program [arg [arg...]]\n\n\
	Runs program with bash-style history written to stdout.\n\
	nnn = number of lines of history to retain.\n\
	If nnn is zero or absent, retain all lines.\n", progname);
	return rc;
	}

 int
main (int argc, char **argv)
{
	char buf[256], *e, *t, *t0;
	int fd[2], i, j, n, n0, p, q;
	static char *signame[16] = { "",
		/*1*/ "Hangup on controlling terminal or death of controlling process",
		/*2*/ "Interrupt from keyboard",
		/*3*/ "Quit from keyboard",
		/*4*/ "Illegal Instruction",
		/*5*/ "Trace/breakpoint trap",
		/*6*/ "Abort signal",
		/*7*/ "Bus error",
		/*8*/ "Floating point exception",
		/*9*/ "Kill -9 signal",
		/*10*/ "Signal 10",
		/*11*/ "Segmentation fault (invalid memory reference)",
		/*12*/ "Signal 12",
		/*13*/ "Broken pipe",
		/*14*/ "Alarmm (SIGALRM)",
		/*15*/ "Termination signal (SIGTERM)"
		};

	n = 0;
	progname = argv[0];
	if (--argc <= 0)
		return usage(1);
	t = *++argv;
	if (*t == '-') {
		if (*++t == '?' && !t[1])
			return usage(0);
		if (!strcmp(t, "-help"))
			return usage(0);
		if (*t >= '0' && *t <= '9') {
			n = (int)strtol(t,&t,0);
			if (*t)
				return usage(1);
			}
		else if (*t != '-' || t[1])
			return usage(1);
		if (!(t = *++argv))
			return usage(1);
		--argc;
		}
	if (pipe(fd)) {
		fprintf(stderr, "%s: pipe failure\n", progname);
		return 2;
		}
	if (!(q = fork())) {
		dup2(fd[0], 0);
		close(fd[0]);
		close(fd[1]);
		execvp(argv[0], argv);
		fprintf(stderr, "Cannot invoke %s\n", argv[0]);
		return 2;
		}
	signal(SIGINT, SIG_IGN);
	signal(SIGPIPE, catch_sigpipe);
	signal(SIGCHLD, catch_sigpipe);
	close(fd[0]);
	p = fd[1];
	using_history();
	if (n)
		stifle_history(n);
	rl_bind_key('\t', rl_insert);	/* treat tab as tab */
	/* history_expansion_char = 0x1b; */	/* escape: treat ! as ! */

	t = t0 = 0;
	while(!pipegone) {
		if (t != t0)
			free(t);
		in_readline = 1;
		if (!(t = readline(0))) {
			pipegone = 1;
			write(p, buf, 0);	/* try to send EOF */
			write(1, "\n", 1);
			break;
			}
		in_readline = 0;
		if (pipegone)
			break;
		n0 = strlen(t);
		if (t[n = n0 - 1] == '\n') {
			if (!n) {
				write(p, t, n0);
				continue;
				}
			t[n] = 0;
			}
		else
			n = -1;
		i = history_expand(t, &e);
		if (i > 0) {
			if (i == 2) {
				n = strlen(e);
				if (n < sizeof(buf)) {
					memcpy(buf, e, n);
					buf[n++] = '\n';
					write(2, buf, n);
					}
				else {
					write(2, e, n);
					if (pipegone)
						break;
					write(2, "\n", 1);
					}
				free(e);
				continue;
				}
			free(t);
			t = e;
			}
		else
			free(e);
		if (!t0 || strcmp(t,t0)) {
			add_history(t);
			if (t0)
				free(t0);
			t0 = t;
			}
		if (i >= 0) {
			n = strlen(t);
			if (n < sizeof(buf)) {
				memcpy(buf, t, n);
				buf[n++] = '\n';
				write(p, buf, n);
				}
			else {
				write(p, t, n);
				if (pipegone)
					break;
				write(p, "\n", 1);
				}
			}
		else {
			if (n >= 0)
				t[n] = '\n';
			write(p, t, n0);
			}
		}
	close(p);
	i = 0;
	do n = wait(&i);
		while(n != -1 && n != q);
	if (j = i & 0xff) {
		if (j < 16)
			fprintf(stderr, "%s\n", signame[j]);
		else
			fprintf(stderr, "Signal %d\n", j);
		return 1;
		}
	return i >> 8;
	}
