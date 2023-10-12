/*
 *	Code adapted from:
 *		The practice of programming/Brian W. Kernighan and Rob Pike,
 *			Addison-Wesley professional computing series, 1999.
 *			ISBN 0-201-61586-X
 *
 */

#ifndef	__LINKED_LIST_VOID_PAYLOAD_HEADER__
#define	__LINKED_LIST_VOID_PAYLOAD_HEADER__

/*
 * define our types
 */
typedef struct LLvNode LLvNode;

struct LLvNode {
	char *key;
	void *value;
	struct LLvNode *next;
};



/* llNewNode: create and initialize data */
LLvNode *llNewNode(char *key, void *value);

/* llPrepend: add newp to front of list */
LLvNode *llPrepend(LLvNode *listp, LLvNode *newp);

/* llAppend: add newp to end of listp */
LLvNode *llAppend(LLvNode *listp, LLvNode *newp);

/* llInsertSort: place the value in the list according to index ordering */
LLvNode *llInsertByIndexOrder(LLvNode *listp, LLvNode *newp);

/* llLookupKey: sequential search for name in listp */
LLvNode *llLookupKey(LLvNode *listp, char *key);

/* llApplyFn: execute fn for each element of listp */
void llApplyFn(LLvNode *listp, void (*fn)(LLvNode*, void*), void *arg);

/* llFree : free all elements of listp */
void llFree(LLvNode *listp, void (*userDeleteFn)(LLvNode*, void*), void *arg);

#endif /*	__LINKED_LIST_VOID_PAYLOAD_HEADER__ */
