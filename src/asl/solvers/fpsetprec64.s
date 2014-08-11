/*  For Solaris (e.g., SunOS 5.11 i386) on an AMD64, with */
/* compilation by cc.  Returns the previous x87 control word. */
	.section	.text,"ax"
	.align 4
	.globl fpsetprec
	.type	 fpsetprec,@function
	.align 16
fpsetprec:
	pushq %rbp
	movq %rsp,%rbp
	subq $16,%rsp
	fstcw -4(%rbp)
	xor %eax,%eax
	movl -4(%rbp),%eax
	movl -4(%rbp),%edx
	andl $-769,%edx
	orl 16(%rbp),%edx
	movl %edx,-4(%rbp)
	fldcw -4(%rbp)
	jmp .L2
	.align 4
.L2:
	leave
	ret
.Lfe1:
	.size	 fpsetprec,.Lfe1-fpsetprec
