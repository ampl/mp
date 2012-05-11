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

#include <stddef.h> /* for size_t */

 typedef struct AVL_Node AVL_Node;
 typedef struct AVL_Tree AVL_Tree;
 typedef struct Element Element;

 typedef int (*AVL_Elcomp)(void*, const Element*, const Element*);
 typedef int (*AVL_Visitor)(void*, Element*);

#ifdef __cplusplus
extern "C" {
#endif

 extern AVL_Tree *AVL_Tree_alloc(void*, AVL_Elcomp, void *(Mallocfunc)(size_t));
 extern AVL_Tree *AVL_Tree_alloc2(void*, AVL_Elcomp, void *(Mallocfunc)(size_t), void (*Free)(void*));
 extern void AVL_Tree_free(AVL_Tree**);
 extern size_t AVL_Tree_size(AVL_Tree*);
 extern Element *AVL_find(const Element *, AVL_Tree*);
 extern Element *AVL_insert(const Element *, AVL_Tree*);
 extern int AVL_visit(void*, AVL_Tree*, AVL_Visitor);
 extern Element *AVL_delete(const Element *, AVL_Tree*);
 extern void *AVL_setv(AVL_Tree *, void*);

 /* The following routines were added March 2011 by David M. Gay (dmg).	*/
 /* AVL_vfind, AVL_vinsert, and AVL_vdelete pass their void* argument	*/
 /* to AVL_Elcomp, rather than the void* value passed to AVL_Tree_alloc,*/
 /* AVL_Tree_alloc2, or AVL_setv.   If the AVL_Node** arguments to	*/
 /* AVL_vfind, AVL_vinsert, AVL_first, and AVL_last are not null, they	*/
 /* are set to point to the AVL_Node associated with the Element found	*/
 /* or, for AVL_vinsert, added.  For AVL_delnode, AVL_next and AVL_prev,*/
 /* the AVL_Node** argument must not be null.   If Free argument to	*/
 /* AVL_Tree_alloc2 is null, it is now taken to be free(). */

 extern Element *AVL_vfind(AVL_Tree*, void*, const Element*, AVL_Node **);
 extern Element *AVL_vinsert(AVL_Tree*, void*, const Element*, AVL_Node **);
 extern Element *AVL_vdelete(AVL_Tree*, void*, const Element*);
 extern Element *AVL_delnode(AVL_Tree*, AVL_Node **);
 extern Element *AVL_first(const AVL_Tree*, AVL_Node**);
 extern Element *AVL_last(const AVL_Tree*, AVL_Node**);
 extern Element *AVL_next(AVL_Node**);
 extern Element *AVL_prev(AVL_Node**);

 /* The following four routines were added November 2011 by dmg to find  */
 /* the first element >= the Element argument or the last element <= the */
 /* Element argument.   The AVL_Node** arguments can be null.  */

 extern Element *AVL_first_ge(const AVL_Tree*, const Element*, AVL_Node**);
 extern Element *AVL_last_le(const AVL_Tree*, const Element*, AVL_Node**);
 extern Element *AVL_vfirst_ge(const AVL_Tree*, void*, const Element*, AVL_Node**);
 extern Element *AVL_vlast_le(const AVL_Tree*, void*, const Element*, AVL_Node**);

 /* The third argument to Avl_Tree_alloc is a malloc-like function that    */
 /* only returns nonzero values.  It should use longjmp to avoid returning */
 /* if no memory is available.  If you are using the AMPL/Solver interface */
 /* library, simply pass mymalloc_ASL for this argument. */

 /* AVL_Tree_alloc returns a pointer to an AVL structure in which it has */
 /* stored its arguments.  The current void* value (possibly a "this"    */
 /* pointer) is passed to the comparision function.   AVL_setv replaces  */
 /* the current such void* value and returns the old value. */

 /* AVL_Tree_size(Tree) returns the number of elements in Tree.		 */

 /* AVL_Visit(Tree,V) calls V once for each element in Tree, from first  */
 /* to last, in the order defined by the AVL_Elcomp associated with Tree,*/
 /* so long as V returns 0.  The return from AVL_visit is the value	 */
 /* returned by the last call on V (or 0 if Tree is empty).		 */

 /* If E is not already in Tree, both AVL_insert(E, Tree) and		 */
 /* AVL_vinsert(E, V, Tree, N) add E to Tree and return 0.  If Tree	 */
 /* contains an Element* E1 that AVL_Elcomp reports to be "equal" to E	 */
 /* (i.e., AVL_Elcomp(v,E,E1) == 0), then AVL_insert() and AVL_vinsert	 */
 /* return E1. */

#ifdef __cplusplus
}
#endif
