/*  For Solaris (e.g., SunOS 5.8 i386) on an x86 (with x87), with */
/* compilation by gcc.  Returns the previous x87 control word. */
	.file	"fpsetprec.s"
	.version	"01.01"
.text
	.align 4
.globl fpsetprec
	.type	 fpsetprec,@function
fpsetprec:
	pushl %ebp
	movl %esp,%ebp
	addl $-4,%esp	/* so -2(%ebp) is above %esp */
	fstcw -2(%ebp)
	xor %eax,%eax
	movw -2(%ebp),%ax
	movw -2(%ebp),%dx
	andw $-769,%dx
	orw 8(%ebp),%dx
	movw %dx,-2(%ebp)
	fldcw -2(%ebp)
	jmp .L2
	.align 4
.L2:
	leave
	ret
.Lfe1:
	.size	 fpsetprec,.Lfe1-fpsetprec
