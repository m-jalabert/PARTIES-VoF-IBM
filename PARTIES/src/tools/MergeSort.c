#include <stdio.h>
#include <stdlib.h>
#include "DataTypes.h"
#include "definitions.h"

/* function prototypes */
MS_STRUCTURE* SortedMerge(MS_STRUCTURE* a, MS_STRUCTURE* b);


void FrontBackSplit(MS_STRUCTURE* source,
          MS_STRUCTURE** frontRef, MS_STRUCTURE** backRef);

/* sorts the linked list by changing next pointers (not MS_DATA) */
void MergeSort(MS_STRUCTURE** headRef)
{
  MS_STRUCTURE* head = *headRef;

  MS_STRUCTURE* a;
  MS_STRUCTURE* b;


  /* Base case -- length 0 or 1 */
  if ((head == NULL) || (head->next == NULL))
  {

    return;
  }
  //printf("Das is nen Fehler sollte auf null zeigen\n");
  //fflush(stdout);

  /* Split head into 'a' and 'b' sublists */
  FrontBackSplit(head, &a, &b);

  /* Recursively sort the sublists */
  MergeSort(&a);
  MergeSort(&b);

  /* answer = merge the two sorted lists together */
  *headRef = SortedMerge(a, b);
}

/* See http://geeksforgeeks.org/?p=3622 for details of this
   function */
MS_STRUCTURE* SortedMerge(MS_STRUCTURE* a, MS_STRUCTURE* b)
{
  MS_STRUCTURE* result = NULL;

  /* Base cases */
  if (a == NULL)
     return(b);
  else if (b==NULL)
     return(a);

  /* Pick either a or b, and recur */
  if (a->MS_DATA <= b->MS_DATA)
  {
     result = a;
     result->next = SortedMerge(a->next, b);
  }
  else
  {
     result = b;
     result->next = SortedMerge(a, b->next);
  }
  return(result);
}

/* UTILITY FUNCTIONS */
/* Split the nodes of the given list into front and back halves,
     and return the two lists using the reference parameters.
     If the length is odd, the extra node should go in the front list.
     Uses the fast/slow pointer strategy.  */
void FrontBackSplit(MS_STRUCTURE* source,
          MS_STRUCTURE** frontRef, MS_STRUCTURE** backRef)
{
  MS_STRUCTURE* fast;
  MS_STRUCTURE* slow;
  if (source==NULL || source->next==NULL)
  {
    /* length < 2 cases */
    *frontRef = source;
    *backRef = NULL;
  }
  else
  {
    slow = source;
    fast = source->next;

    /* Advance 'fast' two nodes, and advance 'slow' one node */
    while (fast != NULL)
    {
      fast = fast->next;
      if (fast != NULL)
      {
        slow = slow->next;
        fast = fast->next;
      }
    }

    /* 'slow' is before the midpoint in the list, so split it in two
      at that point. */
    *frontRef = source;
    *backRef = slow->next;
    slow->next = NULL;
  }
}
