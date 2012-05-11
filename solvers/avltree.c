/****************************************************************
Copyright (C) 2007 David M. Gay

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author disclaims all warranties with regard to this software,
including all implied warranties of merchantability and fitness.
In no event shall the author be liable for any special, indirect or
consequential damages or any damages whatsoever resulting from loss of
use, data or profits, whether in an action of contract, negligence or
other tortious action, arising out of or in connection with the use or
performance of this software.
****************************************************************/

/* AVL_delnode, AVL_first, AVL_last, AVL_next, AVL_prev, AVL_vdelete,
   AVL_vfind, and AVL_vinsert added March 2011 (by dmg). */

#include "avltree.h"
#include <stdlib.h> /* for free() */
#include <string.h>

#ifndef AVL_memgulp
#define AVL_memgulp 256
#endif

 typedef struct AVL_Mblk AVL_Mblk;

 struct
AVL_Node {
	Element *elem;
	AVL_Node *left, *right, *up;
	int height;
	};

 struct
AVL_Mblk {
	AVL_Mblk *next;
	AVL_Node x[AVL_memgulp];
	};

 struct
AVL_Tree {
	AVL_Node *Top;
	AVL_Node *efree;
	AVL_Mblk *mb;
	size_t nelem;
	AVL_Elcomp cmp;
	void *v;
	void *(*Malloc)(size_t);
	void (*Free)(void*);
	};

 AVL_Tree*
AVL_Tree_alloc2(void *v, AVL_Elcomp cmp, void *(*Malloc)(size_t), void (*Free)(void*))
{
	AVL_Mblk *mb;
	AVL_Tree *T;
	AVL_Node *N, *Ne;
	size_t L;

	mb = (AVL_Mblk*)Malloc(L = sizeof(AVL_Tree) + sizeof(AVL_Mblk));
	memset(mb, 0, L);
	T = (AVL_Tree*)(mb + 1);
	T->cmp = cmp;
	T->v = v;
	T->mb = mb;
	T->efree = N = mb->x;
	Ne = N + AVL_memgulp - 1;
	while(N < Ne)
		N = N->left = N + 1;
	N->left = 0;
	T->Malloc = Malloc;
	if (!Free)
		Free = free;
	T->Free = Free;
	return T;
	}

 void
AVL_Tree_free(AVL_Tree **Tp)
{
	AVL_Mblk *mb, *mb1;
	AVL_Tree *T;

	if ((T = *Tp)) {
		*Tp = 0;
		mb1 = T->mb;
		while((mb = mb1)) {
			mb1 = mb->next;
			T->Free(mb);
			}
		}
	}

 size_t
AVL_Tree_size(AVL_Tree *T)
{
	return T->nelem;
	}

 static AVL_Node*
Node_alloc(AVL_Tree *T)
{
	AVL_Mblk *mb;
	AVL_Node *N, *Ne, *Nrv;

	mb = (AVL_Mblk*)T->Malloc(sizeof(AVL_Mblk));
	mb->next = T->mb;
	T->mb = mb;
	N = Nrv = mb->x;
	Ne = N++ + AVL_memgulp - 1;
	T->efree = N;
	while(N < Ne)
		N = N->left = N + 1;
	N->left = 0;
	return Nrv;
	}

 Element *
AVL_vfind(AVL_Tree *T, void *v, const Element *E, AVL_Node **pN)
{
	AVL_Node *N;
	int c;

	for(N = T->Top; N;) {
		if (!(c = (*T->cmp)(v, E, N->elem))) {
			if (pN)
				*pN = N;
			return N->elem;
			}
		if (c < 0)
			N = N->left;
		else
			N = N->right;
		}
	if (pN)
		*pN = 0;
	return 0;
	}

 Element *
AVL_vfirst_ge(const AVL_Tree *T, void *v, const Element *E, AVL_Node **pN)
{
	AVL_Node *N, *Nf;
	int c;

	Nf = 0;
	for(N = T->Top; N;) {
		if (!(c = (*T->cmp)(v, E, N->elem))) {
			Nf = N;
			goto done;
			}
		if (c < 0) {
			Nf = N;
			N = N->left;
			}
		else
			N = N->right;
		}
 done:
	if (pN)
		*pN = Nf;
	if (Nf)
		return Nf->elem;
	return 0;
	}

 Element *
AVL_vlast_le(const AVL_Tree *T, void *v, const Element *E, AVL_Node **pN)
{
	AVL_Node *N, *Nf;
	int c;

	Nf = 0;
	for(N = T->Top; N;) {
		if (!(c = (*T->cmp)(v, E, N->elem))) {
			Nf = N;
			goto done;
			}
		if (c < 0)
			N = N->left;
		else {
			Nf = N;
			N = N->right;
			}
		}
 done:
	if (pN)
		*pN = Nf;
	if (Nf)
		return Nf->elem;
	return 0;
	}

 static void
rebalance(AVL_Tree *T, AVL_Node *N)
{
	AVL_Node *A, *C, *D, *E, *NL, *NR, *NP;
	int d, h, hL, hR;

	for(;; N = NP) {
		hL = hR = 0;
		if ((NL = N->left))
			hL = NL->height + 1;
		if ((NR = N->right))
			hR = NR->height + 1;
		NP = N->up;
		if ((d = hL - hR) > 1) {
			A = NL->left;
			D = NL->right;
			if (!NR) {
				N->left = N->right = 0;
				N->height = 0;
				if (A) {
					if (D) {
						(D->right = N)->up = D;
						NL->height = 2;
						(NL->right = D)->up = NL;
						D->height = 1;
						}
					else {
						(NL->right = N)->up = NL;
						NL->height = 1;
						}
					if (!(NL->up = NP)) {
						T->Top = NL;
						break;
						}
					if (N == NP->left)
						NP->left = NL;
					else
						NP->right = NL;
					}
				else {
					(D->left = NL)->up = D;
					(D->right = N)->up = D;
					NL->height = 0;
					D->height = 1;
					NL->right = 0;
					if (!(D->up = NP)) {
						T->Top = D;
						break;
						}
					if (N == NP->left)
						NP->left = D;
					else
						NP->right = D;
					}
				}
			else if (A->height >= D->height) {
				/* easy rotation */
				(N->left = D)->up = N;
				(NL->right = N)->up = NL;
				NL->height = (N->height = D->height + 1) + 1;
				if (!(NL->up = NP)) {
					T->Top = NL;
					break;
					}
				if (N == NP->left)
					NP->left = NL;
				else
					NP->right = NL;
				}
			else {
				/* hard rotation */
				C = D->left;
				E = D->right;
				if ((NL->right = C))
					C->up = NL;
				--NL->height;
				if ((N->left = E))
					E->up = N;
				(D->left = NL)->up = D;
				(D->right = N)->up = D;
				++D->height;
				N->height = NR->height + 1;
				if (!(D->up = NP)) {
					T->Top = D;
					break;
					}
				if (N == NP->left)
					NP->left = D;
				else
					NP->right = D;
				}
			}
		else if (d < -1) {
			A = NR->right;
			D = NR->left;
			if (!NL) {
				N->left = N->right = 0;
				N->height = 0;
				if (A) {
					if (D) {
						(D->left = N)->up = D;
						NR->height = 2;
						(NR->left = D)->up = NR;
						D->height = 1;
						}
					else {
						(NR->left = N)->up = NR;
						NR->height = 1;
						}
					if (!(NR->up = NP)) {
						T->Top = NR;
						break;
						}
					if (N == NP->left)
						NP->left = NR;
					else
						NP->right = NR;
					}
				else {
					(D->right = NR)->up = D;
					(D->left = N)->up = D;
					NR->height = 0;
					D->height = 1;
					NR->left = 0;
					if (!(D->up = NP)) {
						T->Top = D;
						break;
						}
					if (N == NP->left)
						NP->left = D;
					else
						NP->right = D;
					}
				}
			else if (A->height >= D->height) {
				/* easy rotation */
				(N->right = D)->up = N;
				(NR->left = N)->up = NR;
				NR->height = (N->height = D->height + 1) + 1;
				if (!(NR->up = NP)) {
					T->Top = NR;
					break;
					}
				if (N == NP->left)
					NP->left = NR;
				else
					NP->right = NR;
				}
			else {
				/* hard rotation */
				C = D->right;
				E = D->left;
				if ((NR->left = C))
					C->up = NR;
				--NR->height;
				if ((N->right = E))
					E->up = N;
				(D->right = NR)->up = D;
				(D->left = N)->up = D;
				++D->height;
				N->height = NL->height + 1;
				if (!(D->up = NP)) {
					T->Top = D;
					break;
					}
				if (N == NP->left)
					NP->left = D;
				else
					NP->right = D;
				}
			}
		else {
			if ((h = hL) < hR)
				h = hR;
			if (N->height == h)
				break;
			N->height = h;
			if (!NP) {
				T->Top = N;
				break;
				}
			}
		}
	}

 Element *
AVL_vinsert(AVL_Tree *T, void *v, const Element *E, AVL_Node **pN)
{
	AVL_Node *N, *NP;
	int c;

	if (!(N = T->Top)) {
		if ((N = T->efree))
			T->efree = N->left;
		else
			N = Node_alloc(T);
		N->left = N->right = N->up = 0;
		T->Top = N;
		N->elem = (Element*)E;
		N->height = 0;
		if (pN)
			*pN = N;
		T->nelem = 1;
		return 0;
		}
	for(;;) {
		NP = N;
		if (!(c = (*T->cmp)(v, E, N->elem))) {
			if (pN)
				*pN = N;
			return N->elem;
			}
		if (c < 0) {
			if (!(N = N->left)) {
				if ((N = T->efree))
					T->efree = N->left;
				else
					N = Node_alloc(T);
				N->left = N->right = 0;
				(NP->left = N)->up = NP;
				N->elem = (Element*)E;
				N->height = 0;
				if (!NP->right) {
					rebalance(T,NP);
					}
 ret0:
				if (pN)
					*pN = N;
				++T->nelem;
				return 0;
				}
			}
		else {
			if (!(N = N->right)) {
				if ((N = T->efree))
					T->efree = N->left;
				else
					N = Node_alloc(T);
				N->left = N->right = 0;
				(NP->right = N)->up = NP;
				N->elem = (Element*)E;
				N->height = 0;
				if (!NP->left)
					rebalance(T,NP);
				goto ret0;
				}
			}
		}
	}

 Element*
AVL_delnode(AVL_Tree *T, AVL_Node **pN)
{
	AVL_Node *C, *D, *D1, *N, *NL, *NR, *U;
	Element *E;

	D = N = *pN;
	*pN = 0;
	E = N->elem;

	NR = N->right;
	U = N->up;
	if (!(NL = N->left)) {
		if (!NR) {
			if (!U) {
				T->Top = 0;
				goto done;
				}
			if (N == U->left)
				U->left = 0;
			else
				U->right = 0;
			}
		else {
			if (!U) {
				T->Top = NR;
				NR->up = 0;
				goto done;
				}
			if (N == U->left)
				(U->left = NR)->up = U;
			else
				(U->right = NR)->up = U;
			}
		}
	else if (!NR) {
		if (!U) {
			T->Top = NL;
			NL->up = 0;
			goto done;
			}
		if (N == U->left)
			(U->left = NL)->up = U;
		else
			(U->right = NL)->up = U;
		}
	else {
		if (NL->height <= NR->height) {
			for(D = NL; (D1 = D->right); D = D1);
			C = D->left;
			}
		else {
			for(D = NR; (D1 = D->left); D = D1);
			C = D->right;
			}
		N->elem = D->elem;
		U = D->up;
		if (D == U->left)
			U->left = C;
		else
			U->right = C;
		if (C)
			C->up = U;
		}
	rebalance(T, U);
 done:
	--T->nelem;
	D->left = T->efree;
	T->efree = D;
	return E;
	}

 Element*
AVL_vdelete(AVL_Tree *T, void *v, const Element *E)
{
	AVL_Node *N;
	Element *e;

	if ((e = AVL_vfind(T, v, E, &N)))
		AVL_delnode(T, &N);
	return e;
	}

 Element *
AVL_first(const AVL_Tree * T, AVL_Node **pN)
{
	AVL_Node *N, *N1;

	if (!(N = T->Top))
		return 0;
	while((N1 = N->left))
		N = N1;
	if (pN)
		*pN = N;
	return N->elem;
	}

 Element *
AVL_last(const AVL_Tree *T, AVL_Node **pN)
{
	AVL_Node *N, *N1;

	if (!(N = T->Top))
		return 0;
	while((N1 = N->right))
		N = N1;
	if (pN)
		*pN = N;
	return N->elem;
	}

 Element *
AVL_next(AVL_Node **pN)
{
	AVL_Node *N, *N1;

	N = *pN;
	if ((N1 = N->right)) {
		while((N = N1->left))
			N1 = N;
 ret:
		*pN = N1;
		return N1->elem;
		}
	while((N1 = N->up)) {
		if (N == N1->left)
			goto ret;
		N = N1;
		}
	*pN = 0;
	return 0;
	}

 Element *
AVL_prev(AVL_Node **pN)
{
	AVL_Node *N, *N1;

	N = *pN;
	if ((N1 = N->left)) {
		while((N = N1->right))
			N1 = N;
 ret:
		*pN = N1;
		return N1->elem;
		}
	while((N1 = N->up)) {
		if (N == N1->right)
			goto ret;
		N = N1;
		}
	*pN = 0;
	return 0;
	}

#ifndef AVLTREE_NO_OMIT_NONV

 AVL_Tree*
AVL_Tree_alloc(void *v, AVL_Elcomp cmp, void *(*Malloc)(size_t))
{ return AVL_Tree_alloc2(v, cmp, Malloc, 0); }

 void *
AVL_setv(AVL_Tree *T, void *v)
{
	void *rv = T->v;
	T->v = v;
	return rv;
	}

 static int
avl_visit1(void *v, AVL_Node *N, AVL_Visitor V)
{
	AVL_Node *N1;
	int rv;
 top:
	if ((N1 = N->left))
		avl_visit1(v, N1, V);
	if ((rv = (*V)(v, N->elem)))
		return rv;
	if ((N = N->right))
		goto top;
	return 0;
	}

 int
AVL_visit(void *v, AVL_Tree *T, AVL_Visitor V)
{
	if (!T->Top)
		return 0;
	return avl_visit1(v, T->Top, V);
	}

 Element *
AVL_find(const Element *e, AVL_Tree *T)
{ return AVL_vfind(T, T->v, e, 0); }

 Element *
AVL_insert(const Element *e, AVL_Tree *T)
{ return AVL_vinsert(T, T->v, e, 0); }

  Element *
AVL_delete(const Element *e, AVL_Tree *T)
{ return AVL_vdelete(T, T->v, e); }

  Element *
AVL_first_ge(const AVL_Tree *T, const Element *e, AVL_Node **pN)
{ return AVL_vfirst_ge(T, T->v, e, pN); }

  Element *
AVL_last_le(const AVL_Tree *T, const Element *e, AVL_Node **pN)
{ return AVL_vlast_le(T, T->v, e, pN); }

#endif /* AVLTREE_NO_OMIT_NONV */
