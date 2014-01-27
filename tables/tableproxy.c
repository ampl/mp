/****************************************************************
Copyright (C) 2012 AMPL Optimization, Inc.; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization, Inc. disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
****************************************************************/

/* Proxy table handler */

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <winsock.h>
#include <process.h>
#define socklen_t int
#define SLASH '\\'
#else
#include <unistd.h>
#include <sys/socket.h>
#include <time.h>
#include <netinet/in.h>
#include <netdb.h>
#include <errno.h>
#define SLASH '/'
#endif

#include <stddef.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef STAND_ALONE /*{{*/
#include "asl.h"
#include "avltree.h"
#else /*}{*/
#ifndef _WIN32
#ifdef __hpux
#define NO_SELECT_H
#endif
#ifndef NO_SELECT_H /* for hpux */
#include <sys/select.h>
#endif
#endif
#endif /*}}*/

#include "funcadd.h"
#include "arith.h"	/* for Arith_Kind_ASL and X64_bit_pointers */
#include "tableproxyver.h"

#ifdef X64_bit_pointers
#define Bits "32"
#define Bits1 "64"
#else
#define Bits "64"
#define Bits1 "32"
#endif

/* tableproxy.h { */
#define TABLEPROXY_MAGIC "\nTable_Proxy  \n"

 typedef unsigned int Uint;

 enum { TPH_n_Uint = 10 }; /* number of (adjacent) Uint fields after arith in a TP_Head */

 typedef struct
TP_Head {
	char magic[16];
	Uint arith;
	Uint job;
	Uint TIflags;
	Uint arity;	/* table arity */
	Uint ncols;	/* number of columns */
	Uint nrcols;	/* number of real-valued cols */
	Uint nscols;	/* number of string-valued cols */
			/* Some columns may have both real and string values. */
	Uint nstrings;	/* number of strings in table declaration */
	Uint nrows;	/* For reading, a suggestion for the maximum number rows to return */
			/* at once or the number of rows actually in this batch; */
			/* for writing, the number of rows to write. */
	Uint nrtab;	/* number of values in the real-values table */
	Uint nitab;	/* number of values in the itable (type Uint) */
	Uint tablen;	/* total length of itable, rtable, and stable */

	/* After the header come the real-values table (rtable), the		*/
	/* itable (an array of Uint values), and the string table.  The		*/
	/* rtable has nrcols columns of length nrows, for a total of		*/
	/* nrtab = nrcols * nrows elements.  The itable has nitab elements	*/
	/* and begins with ncols pairs (ir,is), with ir == 0 meaning no		*/
	/* real data for this column ir > 0 meaning real-table column ir-1	*/
	/* is for this column; similarly, is == 0 means no string data for	*/
	/* this column and is > 1 means column is-1 of the itable string	*/
	/* columns is for this column.  Then come nscols * nrows values	for the */
	/* itable string-column array.  The string table starts with the table  */
	/* name, which is followed by nstrings values for the strings in the    */
	/* table declaration, followed by string values for the table.  The     */
	/* elements of the string columns are indexed by elements of the itable */
	/* string-column array, giving offsets in the string table, with offset */
	/* 0 meaning "use the real value" and offset 1 meaning "use the	Missing	*/
	/* data value".	 If nitab > 2*ncols + nscols*nrows, then the itable     */
	/* element following the itable string-column array is the offset of    */
	/* in the string table of "verbose" output for this transaction.	*/

	/* To return an error, we set the TPx_Error bit in job and return	*/
	/* the error message as the one string value (the whole stable).	*/

	} TP_Head;

/* values for TP_Head.job */

enum { TPx_Read = 0, TPx_Write = 1, TPx_Done = 2, TPx_Error = 4,
	TPx_quit = 5, TPx_status = 6 };
/*  } end of tableproxy.h */

 typedef struct ProxProg ProxProg;
 typedef struct MemBlock MemBlock;
 typedef int (*ProxRW)(ProxProg*, void*, Uint);
 typedef const char cchar;

 struct
ProxProg {
	ProxProg *prnext;
	char *name;	/* program or IP addr */
	char *hname;
	ProxRW read, write;
#ifdef _WIN32
	HANDLE hr, hw;
#endif
	int fr, fw;
	int port;	/* -1 for program */
	int proto;
	int remote;
	int rowchunk;
	int verbose;
	};

 struct
MemBlock { union { MemBlock *next; real align; } u; };

 typedef struct
ProxyInfo {
	ProxProg *progs, *remote;
	char *mend, *mnext;
	MemBlock *curmb;
	int proto;
	int remotestarted;
	int verbose;
	} ProxyInfo;

 static ProxyInfo PI0;

#ifndef STAND_ALONE /*{*/
 static void*
Malloc0(AmplExports *ae, size_t L)
{
	void *rv = malloc(L);
	if (!rv) {
		fprintf(Stderr, "malloc(%lu) failure in tableproxy" Bits1 ".dll\n",
			(unsigned long)L);
		exit(1);
		}
	return rv;
	}

 static void*
Malloc1a(AmplExports *ae, size_t L)
{
	MemBlock *mb;
	size_t La;
	void *rv;

	La = PI0.mend - PI0.mnext;
	L = (L + sizeof(real) - 1) & ~(sizeof(real) - 1);
	if (La < L) {
		if (L >= 2048) {
			mb = (MemBlock*)Malloc0(ae, L + sizeof(MemBlock));
			rv = (void*)(mb+1);
			}
		else {
			mb = (MemBlock*)Malloc0(ae, 4096 + sizeof(MemBlock));
			rv = (void*)(mb+1);
			PI0.mnext = (char*)rv + L;
			PI0.mend = (char*)rv + 4096;
			}
		mb->u.next = PI0.curmb;
		PI0.curmb = mb;
		}
	else {
		rv = (void*)PI0.mnext;
		PI0.mnext += L;
		}
	return rv;
	}
#endif /*}*/

 typedef struct IRpair { Uint ir, is; } IRpair;

#ifdef _WIN32 /*{{*/

 static int wsa_started;

 static int
wsa_startup(void)
{
	WORD verreq;
	WSADATA wsaData;
	if (wsa_started)
		return 0;
	verreq = MAKEWORD(2,2);
	if (WSAStartup(verreq, &wsaData))
		return 1;
	wsa_started = 1;
	return 0;
	}

 static int
SocRead(ProxProg *p, void *v, Uint L)
{
	char *s;
	int k, n, rv;

	s = (char*)v;
	k = (int)L;
	for(rv = 0;; s += n) {
		if ((n = recv(p->fr, s, k, 0)) <= 0)
			return -1;
		rv += n;
		if ((k -= n) <= 0)
			break;
		}
	return rv;
	}

 static int
SocWrite(ProxProg *p, void *v, Uint L)
{
	char *s;
	int k, n, rv;

	s = (char*)v;
	k = (int)L;
	for(rv = 0;; s += n) {
		if ((n = send(p->fw, s, k, 0)) <= 0)
			return -1;
		rv += n;
		if ((k -= n) <= 0)
			break;
		}
	return rv;
	}

 static int
Read(ProxProg *p, void *v, Uint L)
{
	DWORD k, n;
	char *s;
	int rv;

	s = (char*)v;
	k = L;
	for(rv = 0;; s += n) {
		if (!ReadFile(p->hr, s, k, &n, 0))
			return -1;
		rv += n;
		if ((k -= n) <= 0)
			break;
		}
	return rv;
	}

 static int
Write(ProxProg *p, void *v, Uint L)
{
	DWORD k, n;
	char *s;
	int rv;

	s = (char*)v;
	k = L;
	for(rv = 0;; s += n) {
		if (!WriteFile(p->hw, s, k, &n, 0))
			return -1;
		rv += n;
		if ((k -= n) <= 0)
			break;
		}
	return rv;
	}
#else /*}{*/

 static int
Read(ProxProg *p, void *v, Uint L)
{
	char *s;
	int k, n, rv;

	s = (char*)v;
	k = (int)L;
	for(rv = 0;; s += n) {
		if ((n = read(p->fr, s, k)) <= 0)
			return -1;
		rv += n;
		if ((k -= n) <= 0)
			break;
		}
	return rv;
	}

 static int
Write(ProxProg *p, void *v, Uint L)
{
	char *s;
	int k, n, rv;

	s = (char*)v;
	k = (int)L;
	for(rv = 0;; s += n) {
		if ((n = write(p->fw, s, k)) <= 0)
			return -1;
		rv += n;
		if ((k -= n) <= 0)
			break;
		}
	return rv;
	}

#define SocRead Read
#define SocWrite Write
#define closesocket close
#endif /*}}*/

 static void
CloseSocket(int sd)
{
	shutdown(sd,2);
	closesocket(sd);
	}

#define TM(len) (*ae->Tempmem)(TI->TMI,len)

 static char*
strcpe(char *t, cchar *f)
{
	while((*t++ = *f))
		++f;
	return t;
	}

 static void
H_init(TP_Head *H, int job, TableInfo *TI)
{
	memset(H, 0, sizeof(TP_Head));
	memcpy(H->magic, TABLEPROXY_MAGIC, sizeof(H->magic));
	H->arith = Arith_Kind_ASL;
	H->job = job;
	if (TI) {
		H->TIflags = TI->flags;
		H->arity = TI->arity;
		H->ncols = TI->ncols;
		}
	}

#ifndef STAND_ALONE /*{{*/

 static ProxProg *
new_PP(AmplExports *ae, ProxyInfo *PI, TableInfo *TI, char *name, int port)
{
	ProxProg *p;
	size_t L;

	L = strlen(name) + 1;
	p = (ProxProg*)Malloc1a(ae, sizeof(ProxProg) + L);
	strcpy(p->name = (char*)(p+1), name);
	if ((p->port = port) > 0) {
		p->prnext = PI->remote;
		PI->remote = p;
		}
	else {
		p->prnext = PI->progs;
		PI->progs = p;
		p->proto = 0;
		}
	return p;
	}

 static void
cleanup(void *v)
{
	MemBlock *mb, *mb1;
	ProxyInfo *Pi;
	ProxProg *p;

	Pi = (ProxyInfo*)v;
	for(p = Pi->progs; p; p = p->prnext) {
		Write(p, "Quit", 5);
#ifdef _WIN32
		CloseHandle(p->hw);
		CloseHandle(p->hr);
#else
		close(p->fw);
		close(p->fr);
#endif
		}
#ifdef _WIN32
	if (Pi->remotestarted)
		WSACleanup();
#endif
	for(mb1 = Pi->curmb; (mb = mb1); ) {
		mb1 = mb->u.next;
		free(mb);
		}
	memset(Pi, 0, sizeof(ProxyInfo));
	}

 static void*
Malloc3(AmplExports *ae, TableInfo *TI, size_t L, cchar *where)
{
	void *rv = malloc(L);
	if (!rv) {
		sprintf(TI->Errmsg = TM(strlen(where) + 64),
			"malloc(%lu) failure in %s.", (unsigned long)L, where);
		}
	return rv;
	}

#ifdef _WIN32 /*{{*/
#define PathName "Path"
 static int
pipe(TableInfo *TI, HANDLE *fd)
{
	SECURITY_ATTRIBUTES S;

	S.nLength = sizeof(S);
	S.lpSecurityDescriptor = 0;
	S.bInheritHandle = FALSE;
	if (CreatePipe(fd,fd+1,&S,4096))
		return 0;
	TI->Errmsg = "CreatePipe failed!";
	return 1;
	}

 static void
inheritable(HANDLE Thisproc, HANDLE *ph)
{
	HANDLE h1 = *ph;

	if (DuplicateHandle(Thisproc, h1, Thisproc, ph, 0, TRUE, DUPLICATE_SAME_ACCESS))
		CloseHandle(h1);
	}

 static char *
make_cmdline(AmplExports *ae, TableInfo *TI, char **a, int sep)
{
	char *e0, *r, *rv;
	char **a0, *s, *s0, *s1;
	int k, k0;

	a0 = a;
	k = 2;
	if (sep)
		while((s = *a++)) {
			s0 = s;
			for(;;) {
				while(*s > ' ') {
					if (*s++ == '"')
						++k;
					}
				if (!*s)
					break;
				k += 2;
				while(*++s <= ' ')
					if (!*s)
						goto break2;
				}
 break2:
			k += s - s0 + 1;
			}
	else
		while((s = *a++))
			k += strlen(s) + 1;
	e0 = 0;
	k0 = k;
	if (!sep) {
		s1 = e0 = GetEnvironmentStrings();
		while(*(s = s1++)) {
			while(*s1++);
			if (s[0] == '=')
				k += s1 - s;
			}
		}
	if (!(r = rv = (char*)Malloc3(ae,TI,k,"make_cmdline")))
		return 0;
	if (k > k0) {
		for(s1 = e0; *(s = s1++); ) {
			while(*s1++);
			if (*s == '=')
				while((*r++ = *s++));
			}
		}
	if (e0)
		FreeEnvironmentStrings(e0);
	if (sep && (s = *a0)) {
		while(*s)
			if (*s++ <= ' ') {
				*r++ = '"';
				for(s = *a0; (*r = *s++); ++r);
				*r++ = '"';
				if (*++a0)
					*r++ = sep;
				break;
				}
		}
	while((s = *a0++)) {
		while((*r = *s++))
			++r;
		*r++ = sep;
		}
	*r = r[1] = 0;
	if (r > rv)
		r[-1] = 0;
	return rv;
	}

#else /*}{*/
#define PathName "PATH"
#endif /*}}*/

 static char *
findprog(AmplExports *ae, char *prog, char *b, size_t buflen)
{
	char *b1, *p, *p0, *s;
	size_t L;
	struct stat sb;
#ifdef _WIN32
#define SEP ';'
	if (!stat(prog,&sb))
		return prog;
	if (*prog && prog[1] == ':')
		goto trystat;
#else /*{*/
#define SEP ':'
	static gid_t gid;
	static uid_t uid;
#ifndef S_IXUSR
#define S_IXUSR 0100
#define S_IXGRP 010
#define S_IXOTH 1
#endif
#ifndef S_IFDIR
#define S_IFDIR 040000
#endif
#endif /*}*/
	for(s = prog; *s; ++s) {
#ifdef _WIN32
		if (*s == '/' || *s == '\\') { /*}*/
 trystat:
#else
		if (*s == '/') {
#endif
			if (stat(prog, &sb))
				return 0;
			return prog;
			}
		}
	if (!(p = getenv(PathName)))
		return p;
#ifndef _WIN32
	if (!gid) {
		gid = getegid();
		uid = geteuid();
		}
#endif
	L = strlen(prog) + 2;
	while(*p) {
		p0 = p;
		while(*p != SEP && *p)
			++p;
		if (p0 == p || (*p0 == '.' && p0+1 == p))
			b1 = prog;
		else if (L + (p-p0) > buflen)
			b1 = 0;
		else {
			for(b1 = b; p0 < p;)
				*b1++ = *p0++;
			while(b1 > b && b1[-1] == SLASH)
				--b1;
			*b1++ = SLASH;
			strcpy(b1, prog);
			b1 = b;
			}
		if (*p == SEP)
			++p;
		if (b1 && !stat(b1, &sb)) {
#ifndef _WIN32
			if (sb.st_mode & S_IFDIR)
				continue;
			if (sb.st_uid == uid) {
				if (!(sb.st_mode & S_IXUSR))
					continue;
				}
			else if (sb.st_gid == gid) {
				if (!(sb.st_mode & S_IXGRP))
					continue;
				}
			else if (!(sb.st_mode & S_IXOTH))
				continue;
#endif
			if (b1 == prog)
				b = prog;
			return b;
			}
		}
	return 0;
	}

 static int
get_sd(AmplExports *ae, TableInfo *TI, char *ipstr, int port, int proto)
{
	int rc, sd;
	size_t L;
	struct hostent *ptrh;
	struct sockaddr_in sab;

	rc = -1;
	memset(&sab, 0, sizeof(sab));
	sab.sin_family = AF_INET;
	sab.sin_port = htons((short)port);
	if (!(ptrh = gethostbyname(ipstr))) {
		L = strlen(ipstr);
		sprintf(TI->Errmsg = TM(L+40), "Bad IP address \"%s\".", ipstr);
		goto ret;
		}
	if ((sd = socket(PF_INET, SOCK_STREAM, proto)) < 0) {
		TI->Errmsg = "socket(PF_INET,...) failed.";
		goto ret;
		}
	memcpy(&sab.sin_addr, ptrh->h_addr, ptrh->h_length);
	if (connect(sd, (struct sockaddr *) &sab, sizeof(sab)) < 0) {
		L = strlen(ipstr);
		sprintf(TI->Errmsg = TM(L+64), "Could not connect to port %d of \"%s\".",
			port, ipstr);
		}
	else
		rc = sd;
 ret:
	return rc;
	}

 static ProxProg*
startremote(AmplExports *ae, TableInfo *TI, ProxyInfo *PI, char *ipstr)
{
	ProxProg *pp;
	char buf[256], *s, *se;
	int port, sd;
	size_t L;
	struct protoent	*ptrp;

	port = 5196;
	L = strlen(ipstr);
	if (L >= sizeof(buf)) {
		sprintf(TI->Errmsg = TM(L + 40), "Oversize 'remote=...' value \"%s\".", ipstr);
		return 0;
		}
	if ((s = strchr(ipstr, '@'))) {
		L = s - ipstr;
		memcpy(buf, ipstr, L);
		buf[L] = 0;
		ipstr = buf;
		port = (int)strtol(++s,&se,10);
		if (port <= 0 || *se) {
			sprintf(TI->Errmsg = TM(strlen(s) + 40), "Bad port number \"%s\".", s);
			return 0;
			}
		}
	if (!PI->remotestarted) {
#ifdef _WIN32
		if (wsa_startup()) {
			TI->Errmsg = "WSAStartup failed.";
			return 0;
			}
#endif
		if (!(ptrp = getprotobyname("tcp"))) {
			TI->Errmsg = "getprotobyname(\"tcp\") failed.";
			return 0;
			}
		PI->proto = ptrp->p_proto;
		PI->remotestarted = 1;
		}
	if ((sd = get_sd(ae, TI, ipstr, port, PI->proto)) < 0)
		return 0;
	pp = new_PP(ae, PI, TI, ipstr, port);
	pp->fr = pp->fw = sd;
	pp->read = SocRead;
	pp->write = SocWrite;
	pp->proto = PI->proto;
	return pp;
	}

 static ProxProg*
startprog(AmplExports *ae, TableInfo *TI, ProxyInfo *PI, char *prog)
{
	ProxProg *pp;
	char *av[3], buf[4096], **ep, *p;
#ifdef _WIN32 /*{{*/
	HANDLE Thisproc, fd0[2], fd1[2], hsave[2];
	PROCESS_INFORMATION Pi;
	STARTUPINFO Si;
	char *cmdline, *d, *env, pbuf[4096], *prog0, *s;
	int b;
	size_t L;
	void(*oldsig)(int);

	/* if prog does not end with ".exe", append ".exe" */
	prog0 = prog;
	d = 0;
	for(s = prog; *s; ++s)
		if (*s == '.')
			d = s;
	if (!d || (s - d) != 4
	 || (d[1] != 'e' && d[1] != 'E')
	 || (d[2] != 'x' && d[2] != 'X')
	 || (d[3] != 'e' && d[3] != 'E')) {
		if ((L = s - prog) >= sizeof(pbuf) - 5)
			return 0;
		memcpy(pbuf, prog, L);
		strcpy(pbuf+L, ".exe");
		prog = pbuf;
		}
#else /*}{*/
	int fd0[2], fd1[2];
	pid_t pid;
#endif /*}}*/

	if (!(p = findprog(ae, prog, buf, sizeof(buf)))) {
		sprintf(TI->Errmsg = TM(strlen(prog) + 32),
			"Could not find \"%s\".", prog);
		return 0;
		}
	av[0] = p;
	av[1] = "--local";
	av[2] = 0;
	ep = ae->ASLdate >= 20111028 && !ae->asl ? (char**)getenv(0) : 0;
#ifdef _WIN32 /*{{*/
	if (!(cmdline = make_cmdline(ae, TI, av, ' ')))
		return 0;
	if (!(env = make_cmdline(ae, TI, ep, 0))) {
 free_cl:
		free(cmdline);
		return 0;
		}
	if (pipe(TI, fd0)) {
 free_e:
		free(env);
		goto free_cl;
		}
	if (pipe(TI, fd1)) {
		CloseHandle(fd0[0]);
		CloseHandle(fd0[1]);
		goto free_e;
		}
	memset(&Si,0,sizeof(Si));
	Si.cb = sizeof(Si);
	Si.wShowWindow = SW_SHOW;
	Thisproc = GetCurrentProcess();
	inheritable(Thisproc, &fd0[1]);
	inheritable(Thisproc, &fd1[0]);
	hsave[0] = GetStdHandle(STD_INPUT_HANDLE);
	hsave[1] = GetStdHandle(STD_OUTPUT_HANDLE);
	SetStdHandle(STD_INPUT_HANDLE, fd1[0]);
	SetStdHandle(STD_OUTPUT_HANDLE, fd0[1]);
	oldsig = signal(SIGINT, SIG_IGN);
	b = CreateProcess(p,cmdline,0,0,1,0,env,0,&Si,&Pi);
	SetStdHandle(STD_INPUT_HANDLE, hsave[0]);
	SetStdHandle(STD_OUTPUT_HANDLE, hsave[1]);
	if (!b) {
		sprintf(TI->Errmsg = TM(strlen(p) + 64),
			"CreateProcess(\"%s\") failure!\nError code %ld.\n",
			p, GetLastError());
		return 0;
		}
	CloseHandle(Pi.hThread);
	CloseHandle(fd1[0]);
	CloseHandle(fd0[1]);
	CloseHandle(Pi.hProcess);
	free(env);
	free(cmdline);
	pp = new_PP(ae, PI, TI, prog0, 0);
	pp->hr = fd0[0];
	pp->hw = fd1[1];
#else /*}{*/
	if (pipe(fd0)) {
 badpipe:
		TI->Errmsg = "pipe failure";
		return 0;
		}
	if (pipe(fd1)) {
		close(fd0[0]);
		close(fd0[1]);
		goto badpipe;
		}
	if ((pid = fork()) == 0) {
		/* child */
		dup2(fd0[1], 1);
		dup2(fd1[0], 0);
		close(fd0[0]);
		close(fd0[1]);
		close(fd1[0]);
		close(fd1[1]);
		signal(SIGINT, SIG_IGN);
		execve(p, av, ep);
		fprintf(Stderr, "execve(\"%s\",...) failed!\n", p);
		exit(1);
		}
	close(fd0[1]);
	close(fd1[0]);
	pp = new_PP(ae, PI, TI, prog, 0);
	pp->fr = fd0[0];
	pp->fw = fd1[1];
#endif /*}}*/
	pp->read = Read;
	pp->write = Write;
	return pp;
	}

 static int
hsetup(AmplExports *ae, TableInfo *TI, ProxyInfo **pPI, ProxProg **pp, int *nsu)
{
	ProxProg *p;
	ProxyInfo *PI;
	char *hname, *ip, *lib, *prog, *s, *s1, *se, **strs, **strs0, **strse, *verb;
	int i, rowchunk, sd, verbose;

	strs = TI->strings;
	strse = strs + TI->nstrings;
	if (strs >= strse || strcmp(strs[0], "tableproxy"))
		return DB_Refuse;
	strs0 = strs;
	ip = lib = prog = verb = 0;
	rowchunk = 512;
	verbose = 0;
	hname = 0;
	if (++strs >= strse) {
		TI->Errmsg = "Two few strings before \":[...]\"; for more on the strings,\n"
		"use the AMPL command\n\n\tprint _handler_desc['tableproxy'];";
		return DB_Error;
		}
	*pPI = PI = (ProxyInfo*)TI->Vinfo;
	verbose = 0;
	for(;;++strs) {
		if (strs >= strse) {
			TI->Errmsg = "No handler name for the proxy table handler.";
			return DB_Error;
			}
		s = s1 = *strs;
		while(*s1 != '=') {
			if (!*s1++)
				goto break2;
			}
		++s1;
		switch(*s) {
		 case 'I':
			if (strncmp(s,"IP=",3)) {
 badv:
				sprintf(TI->Errmsg = TM(strlen(s) + 24),
					"Bad string '%s'.", s);
				return DB_Error;
				}
 ipset:
			if (prog)
				goto notboth;
			ip = s1;
			continue;
		 case 'h':
			if (strncmp(s,"hname=",6))
				goto badv;
			hname = s + 6;
			continue;
		 case 'i':
			if (strncmp(s,"ip=",3))
				goto badv;
			goto ipset;
		 case 'l':
			if (strncmp(s,"lib=",4))
				goto badv;
			if (lib) {
				TI->Errmsg = "Two \"lib=...\" assignments seen; expected at most one.";
				return DB_Error;
				}
			lib = s;
			continue;
		 case 'p':
			if (strncmp(s,"prog=",5))
				goto badv;
			if (ip) {
 notboth:
				TI->Errmsg = "Both \"ip=...\" and \"prog=...\" seen.";
				return DB_Error;
				}
			if (prog) {
				TI->Errmsg = "Two \"prog=...\" assignments seen; expected at most one.";
				return DB_Error;
				}
			prog = s1;
			continue;
		 case 'r':
			if (strncmp(s,"rowchunk=",9))
				goto badv;
			i = (int)strtol(s1,&se,10);
			if (i <= 0 || *se)
				goto badv;
			rowchunk = i;
			continue;
		 case 'v':
			if (strncmp(s,"verbose=",8))
				goto badv;
			verb = s;
			verbose = (int)strtol(s1, &se, 10);
			if (se <= s1)
				goto badv;
			continue;
		 default: goto badv;
		 }
		}
 break2:
	if (ip) {
		for(p = PI->remote; p; p = p->prnext) {
			if (!strcmp(ip, p->name)) {
				if ((sd = get_sd(ae, TI, p->name, p->port, p->proto)) < 0)
					return DB_Error;
				p->fr = p->fw = sd;
				goto found;
				}
			}
		if ((p = startremote(ae, TI, PI, ip)))
			goto found;
		return DB_Error;
		}
	if (!prog)
		prog = "tableproxy" Bits;
	for(p = PI->progs; p; p = p->prnext) {
		if (!strcmp(prog, p->name))
			goto found;
		}
	if (!(p = startprog(ae, TI, PI, prog)))
		return DB_Error;
 found:
	*pp = p;
	if (!hname)
		hname = *strs;
	if (lib)
		*--strs = lib;
	*nsu = strs - strs0;
	p->rowchunk = rowchunk;
	p->verbose = verbose;
	p->hname = hname;
	return 0;
	}

 static size_t
length(char **ps, int n)
{
	int i;
	size_t k;

	k = n;
	for(i = 0; i < n; ++i)
		k += strlen(ps[i]);
	return k;
	}

 static char *
cpall(char *s, char **ps, int n)
{
	int i;
	for(i = 0; i < n; ++i)
		s = strcpe(s, ps[i]);
	return s;
	}

#ifndef ReadGulp
#define ReadGulp 1000
#endif
 enum { Rgulp = ReadGulp };

 static int
Read_ampl_proxy(AmplExports *ae, TableInfo *TI)
{
	DbCol *dbc, *dbc0, *dbce;
	IRpair *irp;
	ProxyInfo *PI;
	ProxProg *p;
	TP_Head *H, H0;
	Uint Lt, i, is, *itab, *itst0, j, nitab, nr, stlen, *z;
	char *s, *stab, **strs, **sv, **sv0, *x[2];
	int k, nc, nc1, nstr, nsu, rc;
	real *rtab, *rtab0;
	size_t L, Ls;

	if ((rc = hsetup(ae,TI,&PI,&p,&nsu)))
		return rc;
	nstr = TI->nstrings - nsu;
	strs = TI->strings + nsu;
	nc = TI->ncols;
	nc1 = nc + TI->arity;
	stlen = strlen(TI->tname) + 1
		+ strlen(p->hname) + 1
		+ length(strs, nstr)
		+ length(TI->colnames, nc1);
	L = sizeof(TP_Head) + stlen;
	rc = DB_Error;
	if (!(H = (TP_Head*)Malloc3(ae,TI,L,"Read_ampl_proxy")))
		goto ret;
	H_init(H, TPx_Read, TI);
	H->nrcols = H->nscols = H->nrtab = H->nitab = 0;
	H->tablen = stlen;
	H->nstrings = nstr;
	H->nrows = p->rowchunk;
	s = stab = (char*)(H+1);
	s = strcpe(s, p->hname);
	s = strcpe(s, TI->tname);
	s = cpall(s, strs, nstr);
	s = cpall(s, TI->colnames, nc1);
	k = p->write(p, H, L);
	free(H);
	if (k < 0) {
		TI->Errmsg = "p->write failure";
		goto ret;
		}
 readmore:
	if ((k = p->read(p, &H0, sizeof(H0))) < 0) {
		TI->Errmsg = "p->read failure";
		goto ret;
		}
	if (H0.job == TPx_Error) {
		if ((Lt = H0.tablen) > 0) {
			TI->Errmsg = (char*)TM(Lt);
			if ((k = p->read(p, TI->Errmsg, Lt)) < 0)
				TI->Errmsg = "p->read_2 error-read failure";
			}
		else
			TI->Errmsg = "Error return with no error message.";
		goto ret;
		}
	if (!(nr = H0.nrows))
		goto done;
	if (H0.ncols != nc) {
		TI->Errmsg = "H0.ncols error";
		goto ret;
		}
	nitab = 2*(H0.ncols + H0.arity) + nr*H0.nscols;
	L = nc1 * sizeof(DbCol);
	Ls = 0;
	if (H0.nscols)
		Ls = (nr * H0.nscols * sizeof(char*) + sizeof(real) - 1) & ~(sizeof(real) - 1);
	if (!(dbc0 = (DbCol*)Malloc3(ae,TI, L + Ls + H0.tablen,"Read_ampl_proxy read alloc")))
		goto ret;
	dbc = dbc0;
	sv0 = (char**)(dbc + nc1);
	rtab = (real*)((char*)sv0 + Ls);
	if (p->read(p, rtab, H0.tablen) < 0) {
		free(dbc0);
		TI->Errmsg = "2nd read in Read_ampl_proxy failed.";
		goto ret;
		}
	itab = (Uint*)(rtab + H0.nrtab);
	stab = (char*)(itab + H0.nitab);
	if (H0.nitab > nitab) {
		s = stab + itab[nitab];
		if (*s)
			printf("%s", s);
		if (!nr) {
			free(dbc0);
			goto done;
			}
		}
	itst0 = itab + 2*nc1 - nr;
	irp = (IRpair*)itab;
	dbce = dbc + nc1;
	memset(dbc, 0, nc1*sizeof(DbCol));
	rtab0 = rtab - nr;
	sv0 -= nr;
	x[0] = 0;
	x[1] = TI->Missing;
	for(; dbc < dbce; ++dbc, ++irp) {
		if ((j = irp->ir))
			dbc->dval = rtab0 + j*nr;
		if ((j = irp->is)) {
			j *= nr;
			dbc->sval = sv = sv0 + j;
			z = itst0 + j;
			for(i = 0; i < nr; ++i)
				sv[i] = (is = z[i]) < 2 ? x[is] : stab + is;
			}
		}
	TI->AddRows(TI, dbc0, nr);
	free(dbc0);
	if (!(H0.job & TPx_Done))
		goto readmore;
 done:
	rc = DB_Done;
 ret:
	if (p->port > 0)
		CloseSocket(p->fr);
	return rc;
	}

 static int
Write_ampl_proxy(AmplExports *ae, TableInfo *TI)
{
	DbCol *db, *db0, *dbe;
	ProxProg *p;
	ProxyInfo *PI;
	TP_Head *H, H0;
	Uint ir, is, *itab, *itc;
	char *Missing, *s, *s0, *s1, **strs, **sv, **sve;
	int nc, nc1, nrc, nsc, nstr, nsu, rc;
	long nr;
	real *r;
	size_t L, itlen, rclen, rtlen, stlen;

	if ((rc = hsetup(ae,TI,&PI,&p,&nsu)))
		return rc;
	rc = DB_Error;
	nc = TI->ncols;
	nc1 = nc + TI->arity;
	db0 = TI->cols;
	dbe = db0 + nc1;
	nstr = TI->nstrings - nsu;
	strs = TI->strings + nsu;
	stlen = strlen(TI->tname) + 1
		+ strlen(p->hname) + 1
		+ length(strs, nstr)
		+ length(TI->colnames, nc1);
	nrc = nsc = 0;
	nr = TI->nrows;
	Missing = TI->Missing;
	for(db = db0; db < dbe; ++db) {
		if (db->dval)
			++nrc;
		if ((sv = db->sval)) {
			++nsc;
			for(sve = sv + nr; sv < sve; ++sv)
				if ((s = *sv) && s != Missing)
					stlen += strlen(s) + 1;
			}
		}
	itlen = 2*nc1 + nsc*nr;
	rtlen = nrc*nr;
	L = sizeof(TP_Head) + rtlen*sizeof(real) + itlen*sizeof(Uint) + stlen;
	if (!(H = (TP_Head*)Malloc3(ae, TI, L, "Write_ampl_proxy")))
		goto ret;
	H_init(H, TPx_Write, TI);
	H->nrcols = nrc;
	H->nscols = nsc;
	H->tablen = L - sizeof(TP_Head);
	H->nstrings = nstr;
	H->nrows = nr;
	H->nrtab = rtlen;
	H->nitab = itlen;
	r = (real*)(H+1);
	itc = (Uint*)(r + rtlen);
	itab = itc + 2*nc1;
	memset(itc, 0, itlen*sizeof(Uint));
	s = s0 = (char*)(itc + itlen);
	s = strcpe(s, p->hname);
	s = strcpe(s, TI->tname);
	s = cpall(s, strs, nstr);
	s = cpall(s, TI->colnames, nc1);
	rclen = nr*sizeof(real);
	ir = is = 0;
	for(db = db0; db < dbe; itc += 2, ++db) {
		if (db->dval) {
			itc[0] = ++ir;
			memcpy(r, db->dval, rclen);
			r += nr;
			}
		if ((sv = db->sval)) {
			itc[1] = ++is;
			for(sve = sv + nr; sv < sve; ++itab, ++sv) {
				if ((s1 = *sv)) {
					if (s1 == Missing)
						*itab = 1;
					else {
						*itab = s - s0;
						s = strcpe(s, s1);
						}
					}
				}
			}
		}
	if (p->write(p, H, L) < 0) {
		TI->Errmsg = "write error in Write_ampl_proxy";
		goto ret;
		}
	if (p->read(p, &H0, sizeof(TP_Head)) < 0) {
		TI->Errmsg = "read error in Write_ampl_proxy";
		goto ret;
		}
	if (memcmp(&H0, TABLEPROXY_MAGIC, sizeof(H->magic))
	 || !(H0.job & (TPx_Done|TPx_Error))) {
		TI->Errmsg = "unexpected reply in Write_ampl_proxy";
		goto ret;
		}
	if (H0.job & TPx_Error) {
		if (p->read(p, TI->Errmsg = TM(H0.tablen), H0.tablen) != H0.tablen)
			TI->Errmsg = "error reading error return in Write_ampl_proxy";
		goto ret;
		}
	if (H0.nitab) {
		if (p->read(p, s = TM(H0.tablen), H0.tablen) != H0.tablen) {
			TI->Errmsg = "error reading verbose output in Write_ampl_proxy";
			goto ret;
			}
		itab = (Uint*)s;
		s = (char*)(itab + H0.nitab) + itab[0];
		if (*s)
			printf("%s", s);
		}
	rc = DB_Done;
 ret:
	if (p->port > 0)
		CloseSocket(p->fr);
	if (H)
		free(H);
	return rc;
	}

#ifdef OTHER_FUNCADD
extern void OTHER_FUNCADD(AmplExports*);
#endif

 void
funcadd(AmplExports *ae)
{
	static char info[] = "tableproxy\n"
	"Proxy table handler for using local " Bits "-bit table handlers\n"
	"and handlers on remote machines (version " TableProxyVersion ").\n"
	"Strings expected before \":[...]\":\n"
	"	'tableproxy' Connection ...\n"
	"where ... are strings for the other handler.  Connection can involve\n"
	"zero or more of\n\n"
	"	'ip=...[@ppp]'\n"
	"	'prog=...'\n"
	"	'hname=...'\n"
	"	'rowchunk=mmm'\n"
	"	'lib=...'\n\n"
	"At most one of 'prog=...' and 'ip=...' may appear.\n"
	"The ... in \"prog=...\" is the desired local program\n"
	"The ... in \"ip=...\" is a remote IP address at which a proxy table handler\n"
	"is running; the ppp in \"@ppp\", if present, is the IP port to use\n"
	"(default = 5196).\n"
	"The ... in 'hname=...' is the handler name (seen in \"display _HANDLERS;\"); if\n"
	"not given, hname is assumed to be the first string in the strings for the other\n"
	"handler (i.e., the ... following Connection).\n"
	"For reading tables, the mmm in 'rowchunk=mmm' is the maximum number of rows for\n"
	"the remote proxy to cache before sending them to the local proxy (default 512).\n"
	"The ... in \"lib=...\" is a shared library in which the remote proxy should\n"
	"look for a suitable handler.  If neither 'ip=' nor 'prog=' appears,\n"
	"'prog=tableproxy" Bits "' is assumed.";
	static int first = 1;

	/* Inform AMPL about the proxy handler. */

	if (first) {
		first = 0;
		at_exit(cleanup, &PI0);
		}
#ifdef OTHER_FUNCADD
	OTHER_FUNCADD(ae);
#endif
	add_table_handler(Read_ampl_proxy, Write_ampl_proxy, info, 0, &PI0);
	}

#else /*} STAND_ALONE {*/


 static char *pfc, *pfc1, *pfce;

 static void
pfc_alloc(size_t *Lp)
{
	char *t;
	size_t L1, L2;
	L1 = pfc1 - pfc;
	if (!pfc)
		L2 = 4096;
	else
		L2 = 2*(pfce - pfc);
	t = (char*)Malloc(L2);
	if (pfc) {
		memcpy(t, pfc, L1);
		free(pfc);
		}
	pfc = t;
	pfc1 = t + L1;
	pfce = t + L2;
	*Lp = pfce - pfc1;
	}

 static int
printf_capture(cchar *fmt, ...)
{
	int n;
	size_t L;
	va_list ap;

	L = pfce - pfc1;
	if (L < 128)
		pfc_alloc(&L);
	va_start(ap, fmt);
	n = vsnprintf(pfc1, L, fmt, ap);
	if (n >= L) {
		do pfc_alloc(&L); while(n >= L);
		n = vsnprintf(pfc1, L, fmt, ap);
		}
	pfc1 += n;
	return n;
	}

 static void*
Malloc1(size_t L)
{
	MemBlock *mb;
	size_t La;
	void *rv;

	La = PI0.mend - PI0.mnext;
	L = (L + sizeof(real) - 1) & ~(sizeof(real) - 1);
	if (La < L) {
		if (L >= 2048) {
			mb = (MemBlock*)Malloc(L + sizeof(MemBlock));
			rv = (void*)(mb+1);
			}
		else {
			mb = (MemBlock*)Malloc(4096 + sizeof(MemBlock));
			rv = (void*)(mb+1);
			PI0.mnext = (char*)rv + L;
			PI0.mend = (char*)rv + 4096;
			}
		mb->u.next = PI0.curmb;
		PI0.curmb = mb;
		}
	else {
		rv = (void*)PI0.mnext;
		PI0.mnext += L;
		}
	return rv;
	}

#ifdef MDEBUG
 static int (*PRintf)(cchar *, ...);
 static int (*FFlush)(FILE*);
 static FILE *STDout;
 int mdbzork, mdbzork1;
 static void
Report(cchar *what, void *v)
{
	PRintf("%d %s: %x\n", ++mdbzork, what, v);
	if (mdbzork == mdbzork1)
		PRintf("");
	FFlush(STDout);
	}
#else
#define Report(a,b) /*nothing*/
#endif

 static char Abits[] = "_" Bits1;

 typedef int (*DbRW)(AmplExports *ae, TableInfo *TI);
 typedef struct
HandlerInfo {
	struct HandlerInfo *nexth;
	DbRW DbRead;
	DbRW DbWrite;
	cchar *info;
	cchar *hname;
	cchar *libname;	 /* name without "_32" or "_64" before final "." */
	cchar *libname_alt; /* name with "_32." or "_64." (NULL if no ".") */
	cchar *libname_found;
	cchar *libname_sought;
	Uint libname_len;
	Uint libname_altlen;
	Uint libname_foundlen;
	Uint libname_soughtlen;
	void *Vinfo;
	int flags;
	} HandlerInfo;

 struct
Element {
	void **v;
	char *vk;
	int elno;
	};

 typedef struct
Tcache {		/* description of a cached column during table read */
	int nr;		/* number of entries (only zero versus nonzero matters) */
	int ns;		/* number of string entries (ditto) */
	real *rc;
	char **sc;
	} Tcache;

 enum {Sblen = 8182};

 typedef struct
Sblock {		/* for string storage during table read */
	struct Sblock *nextsb;
	size_t blen;	/* actual length of sbuf */
	char *sbf;	/* next available byte */
	char *sbe;	/* sbe - sbf == bytes still available in this block */
	/*char sbuf[0];*/
	} Sblock;

 typedef struct
THelp {
	TP_Head tph;
	AVL_Tree *at;		/* for Lookup_w */
	AmplExports *ae;
	HandlerInfo *HI;
	ProxProg *P;
	TMInfo *tmi;
	TableInfo *TI;
	cchar *errmsg;
	char *tabname;
	Tcache *tc;
	Sblock *cursb;
	Sblock *bigsb;
	Sblock *freesb;
	AVL_Tree *sbt;		/* for large Sblocks (for huge strings) */
	Element E;
	long maxrows, nr0;	/* for reading, nr0 = rows so far in this chunk */
				/* for writing, nr0 = number of rows in write table command */
				/* (possibly to be augmented by rows in the existing table) */
	size_t slen;		/* length of strings in this block */
	int maxcols;
	int ndcols, nscols;
	int ntcols;		/* total number of cols: TI->ncols + TI->arity */
	int needswap;
	} THelp;

 typedef union {
	struct sockaddr sa;
	struct sockaddr_in sa4;
#ifndef NO_sockaddr_in6
	struct sockaddr_in6 sa6;
#endif
	} Saddr;

 extern int libload_ASL(AmplExports *ae, cchar *s, int ns, int warn);
 static HandlerInfo *Handlers, **NextHandler = &Handlers;
 static TMInfo TMI;
 static char Loopback[] = "127.0.0.1", Missing[1];
 static char *aalast, *aanext, *curdir;
 static cchar *libname, *libname0, *libname_alt, *libname_found;
 static cchar *libname_sought, *libname0_bits, *libname0_dot, *lnspec;
 static int curdir_len, libname_altlen, libname_len;
 static int libname_foundlen, libname_soughtlen, sigcaught;

 static void
sighand(int sig)
{
	sigcaught = sig;
	}

 static void*
Temp_mem(TMInfo *tmi, size_t L)
{
	TMInfo *tm;

	tm = (TMInfo*)Malloc(sizeof(TMInfo) + L);
	Report("Temp_mem", tm);
	tm->u.prev = tmi->u.prev;
	tmi->u.prev = tm;
	return (void*)(tm+1);
	}

 static void
ignore_addfunc(cchar *name, rfunc f, int type, int nargs, void *funcinfo, AmplExports *ae)
{}

 static void
get_dot(cchar *name, cchar **ndot, cchar **bitsuf)
{
	cchar *dot, *s;

	dot = 0;
	for(s = name; *s; ++s)
		if (*s == '.')
			dot = s;
	*ndot = dot;
	if (dot && dot - name > 3 && dot[-3] == '_'
	 && dot[-2] == Abits[1] && dot[-1] == Abits[2])
		*bitsuf = dot - 3;
	else
		*bitsuf = 0;
	}

#ifdef _WIN32 /*{{*/
 static void
lc_map(char *s)
{
	int c;

	while((c = *s)) {
		if (c >= 'A' && c <= 'Z')
			*s = c + 'a' - 'A';
		++s;
		}
	}

 static int
Abspath(cchar *s)
{
	int c = *s;
	if (((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))
	 && s[1] == ':'
	 && ((c = s[2]) == '\\' || c == '/'))
		return 1;
	return 0;
	}

 static int
has_slash(cchar *s)
{
	int c;

	while((c = *s++))
		if (c == '/' || c == '\\')
			return 1;
	return 0;
	}

#define Getcwd(a,b) GetCurrentDirectory(b,a)

#else /*}{*/

#define lc_map(s) /* do nothing */
#define Abspath(s) (*(s) == '/')
#define Getcwd getcwd

 static int
has_slash(cchar *s)
{
	while(*s)
		if (*s++ == '/')
			return 1;
	return 0;
	}

#endif /*}}*/

 static void
my_add_table_handler(DbRW DbRead, DbRW DbWrite, char *info, int flags, void *Vinfo)
{
	HandlerInfo *H;
	size_t L, L0, L1, Ls;
	char *s;
	cchar *bitsuf, *dot;

	if (!info)
		return;
	while(*info <= ' ')
		if (!*info++)
			return;
	L0 = Ls = 0;
	if (!libname) {
		L0 = libname_len + 1;
		if (libname0_dot)
			L0 += libname_len + (libname0_bits ? -3 : 4);
		Ls = strlen(lnspec);
		L0 += Ls + 1;
		bitsuf = dot = 0;
		get_dot(lnspec, &dot, &bitsuf);
		if (dot)
			L0 += Ls + 4;
		}
	for(L = 0; info[L] > ' '; ++L);
	*NextHandler = H = (HandlerInfo*)Malloc1(sizeof(HandlerInfo) + L0 + L + 1);
	NextHandler = &H->nexth;
	H->nexth = 0;
	H->hname = s = (char*)(H+1);
	memcpy(s, info, L);
	s[L] = 0;
	if (!libname) {
		libname = s += L + 1;
		s = strcpe(s, libname0);
		libname_found = libname;
		libname_foundlen = libname_len;
		libname_alt = 0;
		libname_altlen = 0;
		if (libname0_dot) {
			if (libname0_bits) {
				libname_alt = libname;
				libname = s;
				memcpy(s, libname0, L1 = libname0_bits - libname0);
				s = strcpe(s + L1, libname0_dot);
				libname_altlen = libname_len;
				libname_len -= 3;
				}
			else {
				libname_alt = s;
				memcpy(s, libname0, L1 = libname0_dot - libname0);
				s = strcpe(s + L1, Abits) - 1;
				s = strcpe(s, libname0_dot);
				libname_altlen = libname_len + 3;
				}
			}
		if (!strcmp(libname, libname_sought)) {
			libname_sought = libname;
			libname_soughtlen = libname_len;
			}
		else {
			libname_sought = libname_alt;
			libname_soughtlen = libname_altlen;
			}
		if (Ls) {
			if (bitsuf) {
				libname_alt = s;
				s = strcpe(s, lnspec);
				libname_altlen = Ls;
				L1 = bitsuf - lnspec;
				libname = s;
				memcpy(s, lnspec, L1);
				strcpy(s + L1, dot);
				libname_len = Ls - 3;
				}
			else {
				libname = s;
				s = strcpe(s, lnspec);
				libname_len = Ls;
				libname_altlen = 0;
				if (dot) {
					L1 = dot - lnspec;
					libname_alt = s;
					memcpy(s, lnspec, L1);
					s = strcpe(s+L1, Abits) - 1;
					strcpy(s, dot);
					libname_altlen = Ls + 3;
					}
				}
			}
		}
	while(info[L] && info[L] <= ' ')
		++L;
	H->info = info + L;
	H->DbRead = DbRead;
	H->DbWrite = DbWrite;
	H->flags = flags;
	H->Vinfo = Vinfo;
	H->libname = libname;
	H->libname_len = libname_len;
	H->libname_alt = libname_alt;
	H->libname_altlen = libname_altlen;
	H->libname_sought = libname_sought;
	H->libname_soughtlen = libname_soughtlen;
	H->libname_found = libname_found;
	H->libname_foundlen = libname_foundlen;
	note_libuse_ASL();
	}

 void
auxinfo_ASL(AmplExports *ae)
{
	ae->asl = 0;
	ae->Addfunc = ignore_addfunc;
	ae->Add_table_handler = my_add_table_handler;
	}

 void
af_libnamesave_ASL(AmplExports *ae, cchar *fullname, cchar *name, int nlen)
{
	libname = 0;
	libname0 = fullname;
	libname_len = strlen(fullname);
	get_dot(fullname, &libname0_dot, &libname0_bits);
	}

 static void
swap(void *v, size_t L, Uint n)
{
	char *s, *s1, *s2;
	int t;

	s = (char*)v;
	while(n > 0) {
		--n;
		s1 = s;
		s += L;
		s2 = s - 1;
		do {
			t = *s1;
			*s1++ = *s2;
			*s2-- = t;
			}  while(s1 < s2);
		}
	}

 static int
usage(int rc)
{
	fprintf(rc ? Stderr : stdout, "Usage:\n\n"
	"  %s [--nofork] {job} [port=nnn] [lib [lib...]]\n\n"
	"to start, stop, or inquire about a tableproxy server, using port 5196 or,\n"
	"if given in a port=nnn assignment, port nnn.  If given, {job} must be one of\n\n"
	"\tstart\n"
	"\tstop\n"
	"\trestart\n"
	"\tstatus [ip=IP_addr]\n\n"
	"If a {job} keyword is not given, \"start\" is assumed.  If \"start\" is specified\n"
	"or assumed and lib arguments are given, the specified libraries are loaded\n"
	"initially and table handlers are sought in them; other libraries are loaded\n"
	"on request.\n\n"
	"For \"status\", the optional ip=IP_addr is for checking the status of a remote\n"
	"tableproxy server at the specified (symbolic or numeric) IP address.\n\n"
#ifdef _WIN32
	"--nofork ==> do not detach from the invoking console.\n\n"
#else
	"--nofork ==> neither fork nor ignore SIGHUP, SIGINT, SIGQUIT.\n\n"
#endif
	"When invoked for local use by tableproxy*.dll, the only command-line argument\n"
	"is \"--local\".\n\n", progname);
	return rc;
	}

 static int
Lookup_r(real *r, char **sa, TableInfo *TI)
{ return 0; }

 static int
el_comp(void *v, const Element *a, const Element *b)
{
	char *ak, *ake, *bk, *sa, *sb;
	int j, k, n;
	real d;
	void **va, **vb;
	static int dk[2] = { -1, 1 };

	n = *(int*)v;
	va = a->v;
	vb = b->v;
	ak = a->vk;
	bk = b->vk;
	ake = ak + n;
	while(ak < ake) {
		k = *ak++;
		if (k != *bk++)
			return dk[k];
		if (k) {
			sa = (char*)*va++;
			sb = (char*)*vb++;
			if (sa == sb)
				continue;
			if (sa == Missing)
				return -1;
			else if (sb == Missing)
				return 1;
			if ((j = strcmp(sa, sb)))
				return j;
			}
		else {
			d = *(real*)*va++ - *(real*)*vb++;
			if (d < 0.)
				return -1;
			if (d > 0.)
				return 1;
			}
		}
	return 0;
	}

 static void*
MAalloc0(size_t L)
{
	TMInfo *tm;

	tm = Malloc(L + sizeof(TMInfo));
	Report("MAlloc0", tm);
	tm->u.prev = TMI.u.prev;
	TMI.u.prev = tm;
	return (void*)(tm+1);
	}

 static void*
MAalloc(size_t L)
{
	void *v;

	if (L >= 256)
		return MAalloc0(L);
	L = (L + sizeof(real) - 1) & ~(sizeof(real)-1);
	if (aalast - aanext < L) {
		aanext = (char*) MAalloc0(4096);
		aalast = aanext + 4096;
		}
	v = aanext;
	aanext += L;
	return v;
	}

 static void
MAfree(void* v) {}

 static int
Lookup_w(real *r, char **sa, TableInfo *TI)
{
	AVL_Tree *at;
	DbCol *db, *db0, *dbe;
	Element *e;
	THelp *th;
	char *Ronly, *Sonly, *ek, **sv, *vk;
	int a, j, nd, ns;
	long i, n;
	size_t L, L1, Le;
	void **v;

	th = (THelp*)TI->Private;
	a = TI->arity;
	e = &th->E;
	if (!(at = th->at)) {
		th->at = at = AVL_Tree_alloc2(&th->TI->arity, el_comp, MAalloc, MAfree);
		e->v = (void**)MAalloc(a*(sizeof(void*) + 3));
		e->vk = (char*)(e->v + a);
		e->elno = -1;
		Ronly = e->vk + a;
		Sonly = Ronly + a;
		memset(Ronly, 0, a);
		memset(Sonly, 1, a);
		L = Le = sizeof(Element) + a*sizeof(void*);
		L1 = L + a;
		n = th->nr0;
		db0 = TI->cols;
		dbe = db0 + a;
		for(nd = ns = 0, db = db0; db < dbe; ++db) {
			if (db->dval)
				++nd;
			if (db->sval)
				++ns;
			}
		ek = vk = 0;
		if (!ns)
			ek = vk = Ronly;
		else if (!nd)
			ek = vk = Sonly;
		for(i = 0; i < n; ++i) {
			if (vk == 0) {
				nd = ns = 0;
				for(db = db0; db < dbe; ++db) {
					if ((sv = db->sval) && sv[i])
						++ns;
					else
						++nd;
					}
				Le = L;
				if (!ns)
					ek = Ronly;
				else if (!nd)
					ek = Sonly;
				else {
					ek = 0;
					Le = L1;
					}
				}
			e = (Element*)MAalloc(Le);
			e->v = v = (void**)(e+1);
			e->vk = ek;
			db = db0;
			if (ek == Ronly)
				for(; db < dbe; ++db)
					*v++ = &db->dval[i];
			else if (ek == Sonly)
				for(; db < dbe; ++db)
					*v++ = db->sval[i];
			else {
				e->vk = ek = (char*)(v + a);
				for(; db < dbe; ++db, ++v, ++ek) {
					if ((sv = db->sval) && sv[i]) {
						*ek = 1;
						*v = sv[i];
						}
					else {
						*ek = 0;
						*v = &db->dval[i];
						}
					}
				}
			e->elno = i;
			AVL_insert(e, at);
			}
		e = &th->E;
		}
	ek = e->vk;
	v = e->v;
	j = 0;
	if (r) {
		if (sa)
			for(; j < a; ++j, ++v, ++ek) {
				if ((*v = sa[j]))
					*ek = 1;
				else {
					*ek = 0;
					*v = &r[j];
					}
				}
		else
			for(; j < a; ++j) {
				*ek++ = 0;
				*v++ = &r[j];
				}
		}
	else
		for(; j < a; ++j) {
			*ek++ = 1;
			*v++ = sa[j];
			}
	if ((e = AVL_find(e, at)))
		return e->elno;
	return -1;
	}

 static char *
err_msg(THelp *th, cchar *fmt, ...)
{
	HandlerInfo *H;
	char *s;
	int i;
	size_t L;
	va_list ap;

	va_start(ap, fmt);
	H = th->HI;
	L = strlen(H->hname) + H->libname_soughtlen + strlen(fmt) + 512;
	s = MAalloc(L);
	i = snprintf(s, L, "Table handler %s in %s:\n", H->hname, H->libname_sought);
	vsnprintf(s+i, L-i, fmt, ap);
	return s;
	}

 static int
sb_cmp(void *v, const Element *a, const Element *b)
{
	Sblock *ba, *bb;

	ba = (Sblock*)a;
	bb = (Sblock*)b;
	if (ba->blen < bb->blen)
		return -1;
	if (ba->blen > bb->blen)
		return 1;
	return 0;
	}

 static void
sbcopy(THelp *th, Tcache *tc, long j, char *s)
{
	AVL_Node *N;
	AVL_Tree *sbt;
	Sblock *sb, sb0, *sb1;
	char *t;
	size_t L;

	L = strlen(s) + 1;
	if (L >= (Sblen >> 1)) {
		sb = 0;
		if ((sbt = th->sbt)) {
			sb0.blen = L;
			if ((sb = (Sblock*)AVL_first_ge(sbt, (Element*)&sb0, &N)))
				AVL_delnode(sbt, &N);
			}
		if (!sb) {
			sb = MAalloc0(sizeof(Sblock) + L);
			sb->blen = L;
			}
		sb->nextsb = th->bigsb;
		th->bigsb = sb;
		t = (char*)(sb + 1);
		}
	else {
		if (!(sb = th->cursb) || sb->sbe - sb->sbf < L) {
			if ((sb1 = th->freesb))
				th->freesb = sb1->nextsb;
			else
				sb1 = (Sblock*)MAalloc0(sizeof(Sblock) + Sblen);
			sb1->sbf = (char*)(sb1 + 1);
			sb1->sbe = sb1->sbf + Sblen;
			sb1->nextsb = sb;
			th->cursb = sb = sb1;
			}
		t = sb->sbf;
		sb->sbf += L;
		}
	memcpy(tc->sc[j] = t, s, L);
	th->slen += L;
	}

 static int
rtab_send(THelp *th, TableInfo *TI, int flags)
{
	AVL_Tree *sbt;
	ProxProg *P;
	Sblock *sb, *sb1;
	TP_Head *H;
	Tcache *tc, *tce;
	Uint *cd, *cd0, *it;
	char *s, *s0, **sp, **spe, *t;
	cchar *emsg;
	int nitab, nr, nrtab, ns, ntc, rc;
	long m;
	real *r, *r0;
	size_t L;
	static char verbose_output[] = "\nVerbose table-handler output:\n";

	m = th->nr0;
	th->nr0 = 0;
	if (flags & TPx_Error || (flags & TPx_Done && m == 0)) {
		r0 = 0;
		cd0 = 0;
		L = 0;
		emsg = 0;
		if (flags & TPx_Error) {
			if (!(emsg = th->errmsg))
				emsg = "rtab_send has flags == TPx_Error, but th->errmsg == 0";
			L = strlen(emsg) + 1;
			}
		if (pfc1 > pfc) {
			*pfc1 = 0;
			L += pfc1 - pfc + sizeof(verbose_output) + sizeof(Uint) + 2;
			}
		H = (TP_Head*)Malloc(sizeof(TP_Head) + L);
		H_init(H, flags, TI);
		if (emsg) {
			H->job = TPx_Error;
			s = (char*)(H+1);
			s = strcpe(s, emsg);
			if (pfc1 > pfc) {
				s = strcpe(s-1, verbose_output);
				s = strcpe(s-1, pfc);
				if (s[-2] != '\n') {
					s[-1] = '\n';
					*s = 0;
					}
				pfc1 = pfc;
				}
			goto finish;
			}
		if (pfc1 > pfc) {
			H->nitab = 1;
			it = (Uint*)(H+1);
			s = s0 = (char*)(it + 1);
			}
		goto finish;
		}
	tc = th->tc;
	ntc = th->ntcols;
	tce = tc + ntc;
	for(nr = ns = 0; tc < tce; ++tc) {
		if (tc->nr)
			++nr;
		if (tc->ns)
			++ns;
		}
	nitab = 2*ntc + ns*m;
	if (pfc1 > pfc)
		++nitab;
	nrtab = nr*m;
	L = nitab*sizeof(Uint) + nrtab*sizeof(real) + th->slen;
	if (pfc1 > pfc)
		L += pfc1 - pfc + 1;
	H = Malloc(sizeof(TP_Head) + L);
	H_init(H, TPx_Read | flags, TI);
	H->nrcols = nr;
	H->nscols = ns;
	H->nstrings = 0;
	H->nrows = m;
	H->nrtab = nrtab;
	H->nitab = nitab;
	r = r0 = (real*)(H+1);
	if (nrtab == 0)
		r0 = 0;
	cd = cd0 = (Uint*)(r + nrtab);
	it = cd + 2*ntc;
	s0 = (char*)(cd + nitab);
	s = strcpe(s0, TI->tname);
	nr = ns = 0;
	for(tc = th->tc; tc < tce; ++tc) {
		if (tc->nr) {
			tc->nr = 0;
			memcpy(r, tc->rc, m*sizeof(real));
			r += m;
			*cd++ = ++nr;
			}
		else
			*cd++ = 0;
		if (tc->ns) {
			tc->ns = 0;
			*cd++ = ++ns;
			sp = tc->sc;
			spe = sp + m;
			while(sp < spe) {
				if (!(t = *sp++))
					*it++ = 0;
				else if (t == Missing)
					*it++ = 1;
				else {
					*it++ = s - s0;
					s = strcpe(s, t);
					}
				}
			}
		else
			*cd++ = 0;
		}
 finish:
	if (pfc1 > pfc) {
		*it++ = s - s0;
		*pfc1 = 0;
		s = strcpe(s, pfc);
		pfc1 = pfc;
		}
	H->arith = th->tph.arith;
	H->tablen = L;
	L += sizeof(TP_Head);
	if (th->needswap) {
		if (r0)
			swap(r0, sizeof(real), H->nrtab);
		if (cd0)
			swap(cd0, sizeof(Uint), H->nitab);
		swap(&H->job, sizeof(Uint), TPH_n_Uint);
		}
	P = th->P;
	rc = 0;
	if (P->write(P, H, L) != L) {
		if (TI)
			TI->Errmsg = "write error in rtab_send";
		rc = 1;
		}
	free(H);
	if (rc || flags & (TPx_Done|TPx_Error))
		goto done;
	if ((sb = th->bigsb)) {
		if (!(sbt = th->sbt))
			th->sbt = sbt = AVL_Tree_alloc2(0, sb_cmp, MAalloc, MAfree);
		do AVL_insert((Element*)sb, sbt);
			while((sb = sb->nextsb));
		th->bigsb = 0;
		}
	for(sb = th->cursb; sb; sb = sb1) {
		sb1 = sb->nextsb;
		sb->nextsb = th->freesb;
		th->freesb = sb;
		}
	th->cursb = 0;
 done:
	return rc;
	}

 static void
error_send(THelp *th, TableInfo *TI)
{
	th->errmsg = TI->Errmsg;
	th->nr0 = 0;
	rtab_send(th, TI, TPx_Error);
	}

 static int
add_rows_w(TableInfo *TI, DbCol *cols, long nrows)
{
	/* When writing a table, TI->AddRows should not be called, */
	/* so we give an error return. */
	return 1;
	}

 static int
add_rows_r(TableInfo *TI, DbCol *cols, long nrows)
{
	DbCol *dbc, *dbce;
	Tcache *tc;
	THelp *th;
	char **ps, *s;
	int k, nr, ns;
	long i, j, m;
	real *r;
	size_t L;
	static char msg1[] = "cols[%d].dval == 0, but sval[%ld] == 0";
	static char msg2[] = "cols[%d] has both dval == 0 and sval == 0";

	th = (THelp*)TI->Private;

	/* allocate cache columns if necessary */

	tc = th->tc;
	dbc = cols;
	dbce = dbc + th->ntcols;
	for(nr = ns = 0; dbc < dbce; ++dbc, ++tc) {
		if (dbc->dval && !tc->rc)
			++nr;
		if (dbc->sval && !tc->sc)
			++ns;
		}
	m = th->maxrows;
	if (nr + ns) {
		L = m*(nr*sizeof(real) + ns*sizeof(char*));
		r = (real*)MAalloc(L);
		memset(r, 0, L);
		ps = (char**)(r + nr*m);
		tc = th->tc;
		for(dbc = cols; dbc < dbce; ++dbc, ++tc) {
			if (dbc->dval && !tc->rc) {
				tc->rc = r;
				r += m;
				}
			if (dbc->sval && !tc->sc) {
				tc->sc = ps;
				ps += m;
				}
			}
		}
	j = th->nr0;
	for(i = 0; i < nrows; ++i) {
		tc = th->tc;
		for(dbc = cols; dbc < dbce; ++dbc, ++tc) {
			r = dbc->dval;
			if ((ps = dbc->sval)) {
				if (!(s = ps[i])) {
					if (!r) {
						k = dbc - cols;
						TI->Errmsg = err_msg(th, msg1, k, i);
						return DB_Error;
						}
					goto add_r;
					}
				tc->ns++;
				if (s == Missing) {
					tc->sc[j] = s;
					}
				else
					sbcopy(th, tc, j, s);
				}
			else if (!r) {
				k = dbc - cols;
				TI->Errmsg = err_msg(th, msg2, k);
				return DB_Error;
				}
			else {
 add_r:
				tc->nr++;
				tc->rc[j] = r[i];
				}
			}
		if (++j >= m) {
			th->nr0 = m;
			j = 0;
			if (rtab_send(th, TI, 0))
				return DB_Error;
			}
		}
	th->nr0 = j;
	return 0;
	}

 static long
Adjust_Maxrows(TableInfo *TI, long nmaxrows)
{
	DbCol *dc, *dce;
	THelp *th;
	char **sv, **sv0;
	int nc, nd, ns;
	long nr;
	real *dv, *dv0;

	th = (THelp*)TI->Private;

	nc = TI->ncols + TI->arity;
	dc = TI->cols;
	dce = dc + nc;
	if (!th->maxrows) {
		for(nd = ns = 0; dc < dce; ++dc) {
			if (dc->dval)
				++nd;
			if (dc->sval)
				++ns;
			}
		th->ndcols = nd;
		th->nscols = ns;
		if (nmaxrows < 512)
			nmaxrows = 512;
		dc = TI->cols;
		th->maxrows = nmaxrows;
		}
	else {
		nd = th->ndcols;
		ns = th->nscols;
		if (nmaxrows < (th->maxrows <<= 1))
			nmaxrows = th->maxrows;
		}
	dv = (real*)Temp_mem(th->tmi, nmaxrows*(nd*sizeof(real) + ns*sizeof(char*)));
	sv = (char**)(dv + nmaxrows*nd);
	nr = TI->nrows;
	for(; dc < dce; ++dc) {
		if ((dv0 = dc->dval)) {
			memcpy(dc->dval = dv, dv0, nr*sizeof(real));
			dv += nmaxrows;
			}
		if ((sv0 = dc->sval)) {
			memcpy(dc->sval = sv, sv0, nr*sizeof(char*));
			sv += nmaxrows;
			}
		}
	return TI->maxrows = nmaxrows;
	}

 static char *
s_inc(char **ps)
{
	char *rv, *s;

	rv = s = *ps;
	while(*s++);
	*ps = s;
	return rv;
	}

 static void*
Col_Alloc(TableInfo *TI, int ncol, int sval)
{
	DbCol *db;
	THelp *th;
	long nr;
	size_t L;
	void *v, **vp;

	th = (THelp*)TI->Private;
	if (!(nr = th->maxrows))
		nr = TI->nrows;
	db = TI->cols + ncol;
	if (sval) {
		L = nr*sizeof(char*);
		vp = (void**)&db->sval;
		th->nscols++;
		}
	else {
		L = nr*sizeof(real);
		vp = (void**)&db->dval;
		th->ndcols++;
		}
	memset(*vp = v = (void*)Temp_mem(th->tmi, L), 0, L);
	return v;
	}

 static char *
get_lbuf(cchar *ln, char *lbuf0, size_t lbuf_len, Uint *pls, cchar **nextdir)
{
	Uint ln_len;
	char *lbuf, *s;
	cchar *cdir, *t, *t1;
	int q;
	size_t cdlen, L;

	lbuf = lbuf0;
	ln_len = strlen(ln);
	L = ln_len + 1;
	if (Abspath(ln) || has_slash(ln))
		cdir = 0;
	else {
		if (!(cdir = *nextdir) && (cdir = getenv("ampl_libpath"))) {
			while(*cdir <= ' ' && *cdir)
				++cdir;
			if (!*cdir)
				cdir = 0;
			}
		if (cdir) {
			if ((q = *cdir) == '"' || q == '\'') {
				for(t = ++cdir; *t != q; ++t) {
					if (!*t)
						goto use_curdir;
					}
				if (!(cdlen = t - cdir))
					goto use_curdir;
				while(*++t <= ' ' && *t);
				*nextdir = *t ? t : 0;
				}
			else {
				for(t = cdir;;) {
					if (!*++t) {
						*nextdir = 0;
						break;
						}
					if (*t == '\n') {
						for(t1 = t+1; *t1 && *t1 <= ' '; ++t1);
						*nextdir = *t1 ? t1 : 0;
						break;
						}
					}
				while(t > cdir && t[-1] <= ' ')
					--t;
				cdlen = t - cdir;
				}
			}
		else {
 use_curdir:
			*nextdir = 0;
			if (!curdir && Getcwd(lbuf, lbuf_len)) {
				curdir_len = strlen(lbuf);
				strcpy(curdir = (char*)Malloc(curdir_len + 1), lbuf);
				lc_map(curdir);
				}
			cdir = curdir;
			cdlen = curdir_len;
			}
		}
	if (cdir)
		L += cdlen + 2;
	if (L > lbuf_len)
		lbuf = (char*)Malloc(L);
	s = lbuf;
	if (cdir) {
		memcpy(s, cdir, cdlen);
		s += cdlen;
		*s++ = SLASH;
		}
	strcpy(s, ln);
	lc_map(s);
	*pls = s - lbuf + ln_len;
	return lbuf;
	}

 static int
strcomp(const char *a, const char *b) /* strcmp ignoring case */
{
	int c, d;

	while((c = *a++)) {
		if (!(d = *b++))
			return 1;
		if (c == d)
			continue;
		if (c >= 'A' && c <= 'Z')
			c += 'a' - 'A';
		if (d >= 'A' && d <= 'Z')
			d += 'a' - 'A';
		if (c -= d)
			return c;
		}
	return *b ? -1 : 0;
	}

 static TableInfo *
TI_init(THelp *th, Uint *itab, real *rtab, char *stab)
{
	HandlerInfo *HI, *Hf, **pHI;
	TP_Head *h;
	TableInfo *TI;
	DbCol *dbc;
	Uint i, j, k, len, ls, lseen, nc, nc1, nr, ns, *stc, *sti;
	char *hname, *lbuf, lbuf0[4096], *s, *se, **sp, **st;
	cchar *cannot, *ln, *nextdir, *what;
	int m, ndcols, nscols;
	size_t L;

	th->errmsg = 0;
	th->maxrows = th->nr0 = 0;
	s = stab;
	hname = s_inc(&s);
	th->tabname = s_inc(&s);
	h = &th->tph;
	ns = h->nstrings;
	if (!strncmp(s,"lib=",4)) {
		s += 4;
		ln = s_inc(&s);
		--ns;
		}
	else
		ln = "ampltabl.dll";
	lnspec = ln;
	len = strlen(ln);
	/* look for table hname */
	Hf = 0;
	lseen = 0;
	cannot = nextdir = 0;
	sp = 0; /* silence buggy warning */
	HI = Handlers;
	TI = 0;
	lbuf = lbuf0;
 tryagain:
	for(; HI; HI = HI->nexth) {
		if ((len != HI->libname_len    || strcmp(ln, HI->libname))
		 && (len != HI->libname_altlen || strcmp(ln, HI->libname_alt)))
			continue;
		++lseen;
		if (!strcomp(HI->hname, hname)) {
			Hf = HI;
			if (itab) {
				if (!HI->DbWrite)
					cannot = "write";
				}
			else if (!HI->DbRead)
				cannot = "read";
			break;
			}
		}
	if (!lseen++) {
		pHI = NextHandler;
		lbuf = get_lbuf(ln, lbuf0, sizeof(lbuf0), &ls, &nextdir);
		while((i = libload_ASL(th->ae, libname_sought = lbuf, ls, 2))
				== 1 && nextdir) {
			if (lbuf != lbuf0)
				free(lbuf);
			lbuf = get_lbuf(ln, lbuf0, sizeof(lbuf0), &ls, &nextdir);
			}
		if (!i) {
			if ((HI = *pHI))
				goto tryagain;
			what = "find any table handlers in";
			}
		else
			what = i == 3 ? "find funcadd in" : "load";
		L = ls + 80;
		th->errmsg = se = (char*)Temp_mem(th->tmi, L);
		snprintf(se, L, "Error with proxy table handler:  could not %s \"%s\".",
			what, ln);
		goto ereply;
		}
	if (!(th->HI = Hf) || cannot) {
		L = strlen(hname) + 200;
		if (Hf)
			L += L + Hf->libname_soughtlen;
		th->errmsg = se = (char*)Temp_mem(th->tmi, L);
		m = snprintf(se, L, "Error with proxy table handler: no handler for \"%s\".",
			hname);
		if (Hf)
			m += snprintf(se+m, L-m, "\nDid find %s in \"%s\""
				",\nbut this handler cannot %s tables.\n",
				hname, Hf->libname_sought, cannot);
 ereply:
		th->nr0 = 0;
		rtab_send(th, 0, TPx_Error);
		goto ret;
		}
	nc = h->ncols;
	th->maxcols = th->ntcols = nc1 = nc + h->arity;
	nr = h->nrows;
	L = sizeof(TableInfo) + (nc1+ns)*sizeof(char*) + nc1*sizeof(DbCol);
	if (itab) {
		L +=  nr*h->nscols*sizeof(char*);
		TI = (TableInfo*)Temp_mem(th->tmi, L);
		}
	else {
		L += nc1*(sizeof(real) + sizeof(char*));
		rtab = (real*)Temp_mem(th->tmi, L);
		memset(rtab, 0, nc1*(sizeof(real) + sizeof(char*)));
		sp = (char**)(rtab + nc1);
		TI = (TableInfo*)(sp + nc1);
		}
	th->TI = TI;
	memset(TI, 0, sizeof(TableInfo));
	TI->tname = th->tabname;
	TI->strings = st = (char**)(TI+1);
	for(i = 0; i < ns; ++i)
		*st++ = s_inc(&s);
	TI->colnames = st;
	for(i = 0; i < nc1; ++i)
		*st++ = s_inc(&s);
	TI->Missing = Missing;
	TI->Vinfo = HI->Vinfo;
	TI->TMI = th->tmi;
	TI->nstrings = ns;
	TI->arity = h->arity;
	TI->ncols = nc;
	TI->flags = h->TIflags;
	TI->maxrows = TI->nrows = th->maxrows = nr;
	TI->Private = th;
	TI->AdjustMaxrows = Adjust_Maxrows;
	TI->ColAlloc = Col_Alloc;
	TI->cols = dbc = (DbCol*)(st + nr*h->nscols);
	memset(dbc, 0, nc1*sizeof(DbCol));
	ndcols = nscols = 0;
	if (itab) {
		TI->AddRows = add_rows_w;
		TI->Lookup = Lookup_w;
		stc = itab + 2*nc1;
		for(i = 0; i < nc1; ++i) {
			if ((j = *itab++)) {
				dbc[i].dval = rtab + (j-1)*nr;
				++ndcols;
				}
			if ((j = *itab++)) {
				++nscols;
				--j;
				dbc[i].sval = sp = st + j*nr;
				sti = stc + j*nr;
				for(k = 0; k < nr; ++k) {
					switch(j = sti[k]) {
					  case 0: s = 0; break;
					  case 1: s = Missing; break;
					  default: s = stab + j;
					  }
					sp[k] = s;
					}
				}
			}
		}
	else {
		TI->AddRows = add_rows_r;
		TI->Lookup = Lookup_r;
		for(i = 0; i < nc1; ++i) {
			dbc[i].dval = rtab++;
			dbc[i].sval = sp++;
			}
		}
	th->ndcols = ndcols;
	th->nscols = nscols;
 ret:
	if (lbuf != lbuf0)
		free(lbuf);
	return TI;
	}

 static void
rtab_read(THelp *h, ProxProg *P, char *stab)
{
	HandlerInfo *H;
	TableInfo *TI;
	int k;
	size_t L;

	TI = TI_init(h, 0, 0, stab);
	if (!TI)
		return;
	H = h->HI;
	h->sbt = 0;
	h->cursb  = h->bigsb = h->freesb = 0;
	h->slen = strlen(h->tabname) + 1;
	h->tc = 0;
	L = h->ntcols * sizeof(Tcache);
	memset(h->tc = (Tcache*)MAalloc(L), 0, L);
	switch(k = H->DbRead(h->ae, TI)) {
	 case DB_Done:
		if (rtab_send(h, TI, TPx_Done))
			goto badsend;
		break;
	 case DB_Refuse:
		TI->Errmsg = err_msg(h,
			"\tRefused: missing handler name or database details.");
		goto badsend;
	 default:
		TI->Errmsg = err_msg(h,
			"\tSurprise return %d from table reader.", k);
		/* no break */
	 case DB_Error:
 badsend:
		error_send(h, TI);
		break;
	 }
	}

 static void
wtab_write(THelp *h, ProxProg *P, real *rtab, Uint *itab, char *stab)
{
	HandlerInfo *H;
	TableInfo *TI;
	int rc;

	TI = TI_init(h, itab, rtab, stab);
	if (!TI)
		return;
	h->nr0 = TI->nrows;
	h->at = 0;
	H = h->HI;
	if ((rc = H->DbWrite(h->ae, TI)))
		error_send(h, TI);
	else {
		h->nr0 = 0;
		rtab_send(h, TI, TPx_Write | TPx_Done);
		}
	}

 static int
show_ver(void)
{
	extern char sysdetails_ASL[];
	printf("%s: tableproxy version " TableProxyVersion " (%s, ASL %ld)\n",
		progname, sysdetails_ASL, ASLdate_ASL);
	return 0;
	}

 static void
freeall(THelp *th)
{
	TMInfo *tm, *tmi;

	th->errmsg = 0;
	aanext = aalast = 0;
	tmi = th->tmi;
	tm = tmi->u.prev;
	tmi->u.prev = 0;
	while(tm) {
		Report("freeall", tm);
		tmi = tm->u.prev;
		free(tm);
		tm = tmi;
		}
	}

 static int
quitreply(ProxProg *p, int sd)
{
	TP_Head *H;
	size_t L;
	static char stopmsg[] = "Stopping";

	H = (TP_Head*)Malloc1(L = sizeof(TP_Head) + sizeof(stopmsg));
	H_init(H, TPx_quit, 0);
	strcpy((char*)(H+1), stopmsg);
	H->tablen = sizeof(stopmsg);
	p->write(p, H, L);
	CloseSocket(p->fr);
	return 0;
	}

 static void
show_status(ProxProg *p)
{
	HandlerInfo *HI;
	TP_Head *H;
	char *s, *s0, *se;
	int kh, nh;
	size_t L, L1, Lh;

	nh = 0;
	Lh = 8;
	L = sizeof(TP_Head) + 128;
	for(HI = Handlers; HI; HI = HI->nexth) {
		++nh;
		L1 = strlen(HI->hname);
		if (Lh < L1)
			Lh = L1;
		L += L1 + HI->libname_foundlen;
		}
	H = (TP_Head*)MAalloc(L1 = (nh+1)*(Lh + 2) + L);
	s = s0 = (char*)(H+1);
	se = (char*)H + L;
	kh = Lh;
	if (nh <= 0)
		s += snprintf(s, se-s, "No table handlers loaded.\n");
	else {
		s += snprintf(s, se-s, "%d table handler%s available:\nhandler %*s library\n",
			nh, "s" + (nh == 1), kh-7, "");
		for(HI = Handlers; HI; HI = HI->nexth)
			s += snprintf(s, se-s, "%-*s  %s\n", kh, HI->hname, HI->libname_found);
		}
	H_init(H, TPx_status, 0);
	H->tablen = s - s0;
	p->write(p, H, sizeof(TP_Head) + H->tablen);
	}

 enum { start = 0, stop, status, restart };

 static int
dojob(AmplExports *ae, ProxProg *p, int job, int port, char *s, int sd)
{
	TP_Head H;
	char *ipstr, *s1, *se;
	int j, n, rc;
	size_t L;
	struct hostent *ptrh;
	struct sockaddr_in sab;
	static int jobmap[2] = { TPx_quit, TPx_status };
	static cchar *jobname[2] = { "quit", "status" };

	j = job - 1;
	H_init(&H, jobmap[j], 0);
	ipstr = Loopback;
	p->fr = p->fw = sd;
	rc = 1;
	if (job == status && s) {
		/* look for ip=... */
		if (strncmp(s,"ip=",3) && strncmp(s,"IP=",3)) {
			fprintf(Stderr, "Expected \"ip=...\", not \"%s\"\n", s);
			goto ret;
			}
		ipstr = s += 3;
		for(s1 = s; *s1; ) {
			if (*s1++ == '@') {
				port = (int)strtol(s1, &se, 10);
				if (port <= 0 || *se) {
					fprintf(Stderr, "Bad port value \"%s\" in \"%s\"\n",
						s1, s-3);
					goto ret;
					}
				L = s1 - s;
				ipstr = (char*)Malloc1(L);
				memcpy(ipstr, s, --L);
				ipstr[L] = 0;
				break;
				}
			}
		}
	if (!(ptrh = gethostbyname(ipstr))) {
		fprintf(Stderr, "Bad IP address \"%s\".\n", ipstr);
		goto ret;
		}
	memset(&sab, 0, sizeof(sab));
	sab.sin_family = AF_INET;
	sab.sin_port = htons((short)port);
	memcpy(&sab.sin_addr, ptrh->h_addr, ptrh->h_length);
	if (connect(sd, (struct sockaddr *) &sab, sizeof(sab)) < 0) {
		fprintf(Stderr, "Not running:  could not connect to port %d of \"%s\".\n",
			port, ipstr);
		goto ret;
		}
	if ((n = p->write(p, &H, sizeof(H))) != sizeof(H)) {
		fprintf(Stderr, "write to %s@%d failed.\n", ipstr, port);
		goto ret;
		}
	if ((n = p->read(p, &H, sizeof(H))) != sizeof(H) || H.tablen == 0) {
		fprintf(Stderr, "No reply to %s request.\n", jobname[j]);
		goto ret;
		}
	s1 = (char*)Malloc1(H.tablen+1);
	if ((n = p->read(p, s1, H.tablen)) <= 0) {
		fprintf(Stderr, "Incomplete reply to %s request.\n", jobname[j]);
		goto ret;
		}
	s1[n] = 0;
	printf("%s\n", s1);
	rc = 0;
 ret:
	CloseSocket(p->fr);
	return rc;
	}

 extern int n_badlibs_ASL;
 static char **my_env;

  static char *
my_getenv(cchar *name)
{
	char **e, *s;
	size_t L;

	if (!name)
		return (char*)my_env;
	L = strlen(name);
	for(e = my_env; (s = *e); ++e) {
		if (!strncmp(s,name,L) && s[L] == '=')
			return s + L + 1;
		}
	return 0;
	}

#ifdef _WIN32 /*{*/

 static char*
strcp1(char *t, char *f)
{
	while((*t = *f++))
		++t;
	return t;
	}

 static void
mswinfork(char **argv)
{
	STARTUPINFO Si;
	PROCESS_INFORMATION Pi;
	char *cmdline, *pbuf, pbuf0[2048], *s, *s2;
	int c, i;
	size_t L, La, n;

	L = strlen(progname);
	La = 0;
	for(i = 1; (s = argv[i]); ++i)
		La += strlen(s) + 1;
	pbuf = pbuf0;
	if ((n = La + 2*L + 18) > sizeof(pbuf0)) {
		pbuf = (char*)Malloc(n);
		if (!pbuf) {
			fprintf(Stderr, "%s: Malloc failure in mswinfork()\n", progname);
			exit(1);
			}
		}
	memset(&Si,0,sizeof(Si));
	Si.cb = sizeof(Si);
	Si.wShowWindow = SW_SHOWMINIMIZED;
	for(s = pbuf, s2 = progname; (c = *s2++); s++) {
		if (c == '/')
			c = '\\';
		*s = c;
		}
	*s = 0;
	if (s - pbuf > 4 && s[-4] != '.') {
		strcpy(pbuf+L, ".exe");
		L += 4;
		}
	cmdline = pbuf + ++L;
	cmdline[0] = '"';
	s = strcp1(cmdline+1, pbuf);
	s = strcp1(s, "\" --nofork");
	for(i = 1; (s2 = argv[i]); ++i) {
		*s++ = ' ';
		s = strcp1(s, s2);
		}
	c = CreateProcess(pbuf,cmdline,0,0,1,0,0,0,&Si,&Pi);
	if (!c) {
		fprintf(Stderr, "%s: CreateProcess(\"%s\") failed in mswinfork()\n",
			progname, pbuf);
		exit(1);
		}
	CloseHandle(Pi.hThread);
	CloseHandle(Pi.hProcess);
	exit(0);
	}
#endif /*}*/

 int
main(int argc, char **argv, char **arge)
{
	ASL *asl;
	ProxProg P;
	Saddr sab, sac;
	THelp th;
	Uint *itab;
	char **av, *s, *se, *stab, *t;
	fd_set fds;
	ProxRW Rd, Wr;
	int dofork, i, job, n, port, proto, rc, rsdf, rstrt, sd, sdi, sdmx;
	real *rtab;
	size_t L;
	socklen_t soclen;
	struct protoent	*ptrp;
	void (*sighandler)(int);

#ifdef _WIN32 /*{{*/
	char **av0 = argv, **prs = 0;
#define NO_REUSEADDR
#else /*}{*/
#ifndef Use_sleep
	struct timespec wrem, wtime;
#endif
#endif /*}}*/

	/*Sleep(40000);*/	/*DEBUG*/
	progname = *argv;
	my_env = arge;
	port = 5196;
	dofork = rc = rsdf = 1;
	Stderr_init_ASL();
	proto = rstrt = 0;
	sighandler = sighand;
 nextopt:
	if ((s = *++argv) && *s == '-')
		switch(*++s) {
		 case '?':
		 case 'h':
			return usage(s[1] != 0);
		 case '-':
			if (!strcmp(++s, "local")) {
#ifdef DEBUG_FIFO
				char dbuf[64];
				int dbf, dblen;

				dblen = snprintf(dbuf, sizeof(dbuf),
					"echo to debug.fifo to continue /proc/%lu\n",
					(unsigned long)getpid());
				write(2, dbuf, dblen);
				if ((dbf= open("debug.fifo", O_RDONLY)) >= 0) {
					read(dbf, dbuf, sizeof(dbuf));
					close(dbf);
					}
#endif
				if (argv[1])
					goto baduse;
				port = 0;
				s = 0;
				break;
				}
			if (!strcmp(s, "help"))
				return usage(0);
			if (!strcmp(s,"version"))
				return show_ver();
			if (!strcmp(s,"nofork")) {
				dofork = 0;
				goto nextopt;
				}
			if (!*s) {
				s = *++argv;
				break;
				}
			return usage(1);
		 case 'v':
			if (!s[1])
				return show_ver();
		 default:
 baduse:
			return usage(1);
		 }
	job = start;
	if (s) {
		if (s[0] == 's' && s[1] == 't') {
			if (!strcmp(s+2,"art"))
				goto nextarg;
			else if (!strcmp(s+2, "op")) {
				job = stop;
				goto nextarg;
				}
			else if (!strcmp(s+2, "atus")) {
				job = status;
 nextarg:
				s = *++argv;
				}
			}
		else if (!strcmp(s, "restart")) {
			job = rstrt = restart;
#ifdef _WIN32
			prs = argv;
#endif
			goto nextarg;
			}
		}
	if (s) {
		if (!strncmp(s,"port=",5)) {
			port = (int)strtol(s += 5, &se, 10);
			if (port <= 0 || *se) {
				fprintf(Stderr, "Bad \"port=\" value \"%s\"\n", s);
				return usage(1);
				}
			s = *++argv;
			}
		}
	i_option_ASL = "";	/* do not load any table handlers */
	if (!job && s) {
		L = 0;
		for(av = argv; s; s = *++argv)
			L += strlen(s) + 1;
		if (L) {
			i_option_ASL = t = (char*)Malloc(L);
			while((s = *av++)) {
				t = strcpe(t, s);
				t[-1] = '\n';
				}
			t[-1] = 0;
			}
		}
	sd = sdi = 0;
	memset(&P, 0, sizeof(P));
	if (port) {
#ifdef _WIN32
		if (dofork && !job)
			mswinfork(av0);
		if (wsa_startup()) {
			fprintf(Stderr, "WSAStartup failed.\n");
			return 1;
			}
#endif
 more_restart:
		Rd = SocRead;
		Wr = SocWrite;
		if (!(ptrp = getprotobyname("tcp"))) {
			fprintf(Stderr, "Can't map TCP to protocol number.");
			return 1;
			}
		proto = ptrp->p_proto;
		if ((sd = socket(PF_INET, SOCK_STREAM, proto)) < 0) {
			fprintf(Stderr, "socket(...) failed.\n");
			return 1;
			}
		sdmx = sd + 1;
		if (job) {
			rsdf = dofork;
			dofork = 0;
			goto job1;
			}
#ifndef NO_REUSEADDR
		i = 1;
		if (setsockopt(sd, SOL_SOCKET, SO_REUSEADDR, &i, sizeof(i)))
			fprintf(Stderr, "setsockopt(...) failed.\n");
#endif
		memset(&sab, 0, sizeof(sab));
		sab.sa4.sin_family = AF_INET;
		sab.sa4.sin_port = htons((short)port);
		if (bind(sd, (struct sockaddr*)&sab.sa4, sizeof(struct sockaddr))) {
			perror("bind(...) failed");
			goto done;
			}
		if (listen(sd, 128)) {
			fprintf(Stderr, "listen(...) failed\n");
			goto done;
			}
		FD_ZERO(&fds);
		if (rstrt)
			goto more_restart2;
		}
	else {
		dofork = 0;
		sdmx = 2;
		Rd = Read;
		Wr = Write;
		P.fw = 1;
#ifdef _WIN32
		P.hr = GetStdHandle(STD_INPUT_HANDLE);
		P.hw = GetStdHandle(STD_OUTPUT_HANDLE);
#endif
		}
 job1:
	P.read = Rd;
	P.write = Wr;
	th.P = &P;
	asl = (ASL*)Malloc(sizeof(ASL));
	memset(asl, 0, sizeof(ASL));
	asl->p.need_funcadd_ = 1;
	/* get initial table handlers; we may load more later */
	func_add(asl);
	if (n_badlibs_ASL)
		goto done;
	th.ae = asl->i.ae;
	th.ae->Getenv = my_getenv;
	th.ae->Tempmem = Temp_mem;
#ifdef MDEBUG
	PRintf = printf;
	FFlush = fflush;
	STDout = stdout;
#endif
	th.ae->PrintF = printf_capture;
 more_restart2:
#ifndef _WIN32 /*{*/
	signal(SIGCHLD, SIG_IGN);
	if (dofork) {
		sighandler = SIG_IGN;
		if (fork())
			_exit(0);
		setsid();
		if (fork())
			_exit(0);
		freopen("/dev/null", "r", stdin);
		freopen("/dev/null", "w", stdout);
		freopen("/dev/null", "w", stderr);
		}
	signal(SIGHUP,  sighandler);
	signal(SIGQUIT, sighandler);
	signal(SIGTERM, sighand);
#endif /*}*/
	signal(SIGINT, sighandler);
	th.tmi = &TMI;
	th.at = 0;
	if (job) {
		if (rstrt) {
			dojob(asl->i.ae, &P, stop, port, s, sd);
			CloseSocket(sd);
			dofork = rsdf;
			job = start;
#ifdef _WIN32 /*{{*/
			Sleep(2000);	/* 2 seconds */
#else /*}{*/
#ifdef Use_sleep
			sleep(2);
#else
			wtime.tv_sec = 2;
			wtime.tv_nsec = 0;
			nanosleep(&wtime, &wrem);
#endif
#endif /*}}*/
			printf("Restarting\n");
#ifdef _WIN32
			if (dofork) {
				*prs = "start";
				mswinfork(av0);
				}
#endif
			goto more_restart;
			}
		rc = dojob(asl->i.ae, &P, job, port, s, sd);
		goto done;
		}
	for(;;) {
		if (sigcaught)
			break;
		if (port) {
			FD_SET(sd, &fds);
			n = select(sdmx, &fds, 0, 0, 0);
			if (sigcaught)
				break;
			if (n <= 0 || !FD_ISSET(sd, &fds))
				continue;
			soclen = sizeof(sac);
			sdi = accept(sd, &sac.sa, &soclen);
			if (sdi < 0) {
				perror("error in accept()");
				goto done;
				}
			P.fr = P.fw = sdi;
			}
		n = Rd(&P, &th.tph, sizeof(TP_Head));
		if (n < sizeof(TP_Head))
			break;
		if (strcmp(th.tph.magic, TABLEPROXY_MAGIC))
			break;
		if ((th.needswap = th.tph.arith != Arith_Kind_ASL))
			swap(&th.tph.job, sizeof(Uint), TPH_n_Uint);
		switch(th.tph.job) {
		  case TPx_Read:
		  case TPx_Write:
			break;
		  case TPx_quit:
			rc = quitreply(&P, sdi);
			goto done;
		  case TPx_status:
			show_status(&P);
			/* no break */
		  default:
			goto loopend1;
		  }
		rtab = (real*)Temp_mem(th.tmi, th.tph.tablen);
		itab = (Uint*)(rtab + th.tph.nrtab);
		stab = (char*)(itab + th.tph.nitab);
		n = th.tph.tablen;
		s = (char*)rtab;
		for(;;) {
			i = Rd(&P, s, n);
			if (i <= 0)
				goto loopend;
			if ((n -= i) <= 0)
				break;
			s += i;
			}
		if (th.needswap && stab > (char*)rtab) {
			swap(rtab, sizeof(real), th.tph.nrtab);
			swap(itab, sizeof(Uint), th.tph.nitab);
			}
		if (th.tph.job == TPx_Read)
			rtab_read(&th, &P, stab);
		else
			wtab_write(&th, &P, rtab, itab, stab);
 loopend:
		freeall(&th);
 loopend1:
		if (port)
			CloseSocket(sdi);
		}
 done:
#ifdef _WIN32
	if (port)
		WSACleanup();
#endif
	return rc;
	}

#endif /*}} STAND_ALONE*/
