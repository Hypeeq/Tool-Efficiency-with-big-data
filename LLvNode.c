/*
 *	Code adapted from:
 *		The practice of programming/Brian W. Kernighan and Rob Pike,
 *			Addison-Wesley professional computing series, 1999.
 *			ISBN 0-201-61586-X
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "LLvNode.h"


/*
 * llNewNode: create and initialize data
 */
LLvNode *
llNewNode(char *key, void *value)
{
	LLvNode *newp;

	newp = (LLvNode *) malloc(sizeof(LLvNode));

	/* assign data within new node */
	newp->key = key;
	newp->value = value;

	/* make sure we point at nothing */
	newp->next = NULL;

	return newp;
}


/*
 * llAppend: add newp to end of listp
 *
 * as above, we always return the value that
 * should be the new head of the list 
 */
LLvNode *
llAppend(LLvNode *listp, LLvNode *newp)
{
	LLvNode *p;

	if (listp == NULL)
		return newp;

	for (p = listp; p->next; p = p->next)
		;

	p->next = newp;
	return listp;
}


/*
 * llPrepend: add newp to front of list
 *
 * by returning ourselves, we give back the
 * value that should be the new head of the list
 */
LLvNode *
llPrepend(LLvNode *listp, LLvNode *newp)
{
	newp->next = listp;
	return newp;
}


/* llLookupKey: sequential search for key in listp */
LLvNode *
llLookupKey(LLvNode *listp, char *key)
{
	for ( ; listp != NULL; listp = listp->next) {
		if (strcmp(key, listp->key) == 0) {
			return listp;
        }
    }

	return NULL; /* no match */
}


/* llApplyFn: execute fn for each element of listp */
void
llApplyFn(LLvNode *listp, void (*fn)(LLvNode*, void*), void *arg)
{
	for ( ; listp != NULL; listp = listp->next)
		(*fn)(listp, arg); /* call the function */
}


/* llFree : free all elements of listp */
void
llFree(LLvNode *listp, void (*userDeleteFn)(LLvNode*, void*), void *arg)
{
	LLvNode *next;

	for ( ; listp != NULL; listp = next) {

		/** hang on to the next pointer */
		next = listp->next;

		/** call the user function, if provided */
		if (userDeleteFn != NULL)
			(*userDeleteFn)(listp, arg);

		/** free the list node itself */
		free(listp);
	}
}

