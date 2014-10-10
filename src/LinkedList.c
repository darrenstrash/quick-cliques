/* 
    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version. 
 
    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
    GNU General Public License for more details. 
 
    You should have received a copy of the GNU General Public License 
    along with this program.  If not, see <http://www.gnu.org/licenses/> 
*/

#include"LinkedList.h"
#include"misc.h"
#include"MemoryManager.h"
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>

/*! \file LinkedList.c

    \brief a doubly-linked list data structure, with location-aware
           functionality.

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly

    This file contains two types of functions, those that operate
    on Links and those that operate on LinkedLists. These have not
    been optimized for speed or memory usage. 

    Sentinels are used as the head and tail to mark the first Link
    and the last Link in the LinkedList. This makes the logic for
    adding and removing Links more succinct.

*/

/** \name Functions on Link structures
    \brief These functions operate on Link structures.
 */

//@{

/*! \brief tells if a given Link is the head sentinel.

    \param list A link structure.

    \return true if the link is the head sentinel, false otherwise.
*/

int isHead(Link* list)
{
#ifdef DEBUG
    printf("isHead...\n");
#endif

    assert(list != NULL);

    return (list->prev == NULL);
}

/*! \brief tells if a given Link is the tail sentinel.

    \param list A link.

    \return true if the link is the tail sentinel, false otherwise.
*/

int isTail(Link* list)
{
#ifdef DEBUG
    printf("isTail...\n");
#endif

    assert(list != NULL);

    return (list->next == NULL);
}

/*! \brief location-aware function to add a link after a given link.

    \param list The link that we want to add our data after.

    \param data A peice of data to put in the added link.

    \return A pointer to the Link that was added after list.
*/

Link* addAfter(Link* list, void* data)
{
#ifdef DEBUG
    printf("addAfter...\n");
#endif

    assert(list != NULL);
    assert(list->next != NULL);

    Link* newLink = Malloc(sizeof(Link));

    newLink->data = data;

    newLink->next = list->next;
    newLink->prev = list;

    list->next->prev = newLink;
    list->next = newLink;

#ifdef SMARTLENGTH
    newLink->linkedList = list->linkedList;
    newLink->linkedList->length++;
#endif

    return newLink;
}

/*! \brief location-aware function to add a link before a given link.

    \param list The link that we want to add our data before. 

    \param data A peice of data to put in the added link.

    \return A pointer to the Link that was added after list.
*/

Link* addBefore(Link* list, void* data)
{
#ifdef DEBUG
    printf("addBefore...\n");
#endif

    assert(list != NULL);
    assert(list->prev != NULL);

    Link* newLink = Malloc(sizeof(Link));

    newLink->data = data;

    newLink->next = list;
    newLink->prev = list->prev;

    list->prev->next = newLink;
    list->prev = newLink;

#ifdef SMARTLENGTH
    newLink->linkedList = list->linkedList;
    newLink->linkedList->length++;
#endif

    return newLink;
}

/*! \brief delete the given link, and return its data. 

    \param list The link that we want to get rid of.

    \return the data that was in the link, if it was 
            allocated by you, you need to free it..
*/

void* delete(Link* list)
{
#ifdef DEBUG
    printf("delete...\n");
#endif
    assert(list != NULL);
    assert(list->next != NULL);
    assert(list->prev != NULL);

    void* data = list->data;

    Link* linkToFree = removeLink(list);

    Free(linkToFree);

    return data;
}

/*! \brief location-aware method to add a link before another link. 

    \param list The link that we want to add a link before.

    \param newLink The Link to be added after list.
*/

void addLinkBefore(Link* list, Link* newLink)
{
    assert(list != NULL);
    assert(list->prev != NULL);
    assert(newLink != NULL);

    newLink->next = list;
    newLink->prev = list->prev;

    newLink->next->prev = newLink;
    newLink->prev->next = list->prev;

#ifdef SMARTLENGTH
    newLink->linkedList = list->linkedList;
    newLink->linkedList->length++;
#endif
    
}

/*! \brief location-aware method to remove a link, and return it.

    \param list The link that we want remove and return.

    \return the removed link
*/

Link* removeLink(Link* list)
{
#ifdef DEBUG
    printf("removeLink...\n");
#endif
    assert(list != NULL);
    assert(list->next != NULL);
    assert(list->prev != NULL);

    list->next->prev = list->prev;
    list->prev->next = list->next;

    list->next = NULL;
    list->prev = NULL;

#ifdef SMARTLENGTH
    list->linkedList->length--;
    list->linkedList = NULL;
#endif

    return list;
}

//@}

/** \name Functions on LinkedList structures
    \brief These functions operate on LinkedList structures.
 */

//@{

/*! \brief create a new empty linked list

    \return the created linked list
*/

LinkedList* createLinkedList(void)
{
    LinkedList* linkedList = Malloc(sizeof(LinkedList));

    linkedList->head = Malloc(sizeof(Link));
    linkedList->tail = Malloc(sizeof(Link));

    linkedList->head->prev = NULL;
    linkedList->head->next = linkedList->tail;
    linkedList->head->data = (void*) 0xDEAD0000;

    linkedList->tail->prev = linkedList->head;
    linkedList->tail->next = NULL;
    linkedList->tail->data = (void*) 0xDEADFFFF;
#ifdef SMARTLENGTH
    linkedList->length = 0;
    linkedList->head->linkedList = linkedList;
    linkedList->tail->linkedList = linkedList;
#endif

    return linkedList;
}

/*! \brief destroy a linked list

    If you allocated data that is in each link, then
    this will cause a memory leak for you.

    \see destroyLinkedListWithClean

    \param linkedList The linked list to destroy.
*/

void destroyLinkedList(LinkedList* linkedList)
{
    Link* curr = linkedList->head;

    while(curr != NULL)
    {
        Link* currNext = curr->next;
        Free(curr);
        curr = currNext;
    }

    Free(linkedList);
}

/*! \brief destroy a linked list and run a clean function
           on the data in each link.

    \param linkedList The linked list to destroy.

    \param clean A pointer to a function that cleans the data in the links.
*/

void destroyLinkedListWithClean(LinkedList* linkedList, void (*clean)(void*))
{
    Link* curr = linkedList->head;

    while(curr != NULL)
    {
        Link* currNext = curr->next;
        clean(curr->data);
        Free(curr);
        curr = currNext;
    }

    Free(linkedList);
}

/*! \brief copy a linked list

    \param destination copy the linked list here

    \param source copy this linked list
*/

void copyLinkedList(LinkedList* destination, 
                    LinkedList* source)
{
    assert(destination != NULL && source != NULL);

    Link* curr = source->head->next;

    while(!isTail(curr))
    {
        addLast(destination, curr->data);
        curr = curr->next;
    }
}

/*! \brief Compare two linked lists to see if they are equal.

    \param list1 A linked list.

    \param list2 Another linked list.

    \param comparator A function to compare data in the links copy this linked list.

    \return true if the input linked lists have the same data in the same order.
*/

int equal( LinkedList* list1, 
           LinkedList* list2, 
           int (*comparator)(void*,void*))
{
    assert(list1 != NULL && list2 !=NULL);

    Link* curr1 = list1->head->next;
    Link* curr2 = list2->head->next;

    while(!isTail(curr1) && !isTail(curr2))
    {
        if(comparator(curr1->data, curr2->data) == 0)
        {
            curr1 = curr1->next;
            curr2 = curr2->next;
        }
        else if(comparator(curr1->data, curr2->data) > 0)
        {
            return 0;
            curr2 = curr2->next;
        }
        else
        {
            return 0;
            curr1 = curr1->next;
        }
    }

    return (isTail(curr1) && isTail(curr2));
}

/*! \brief decide if a linked list contains a piece of data.

    \param linkedList A linked list.

    \param data The data that we want to look for in linkedList.

    \param comparator A function that returns 0 when two data elements are equal.

    \return true if linkedList contains data.
*/

int contains(LinkedList* linkedList, void* data, int (*comparator)(void*,void*))
{
    assert(linkedList != NULL);

    Link* curr = linkedList->head->next;

    while(!isTail(curr))
    {
        if(comparator(curr->data, data) == 0)
        {
            return 1;
        }
        curr = curr->next;
    }

    return 0;
}

/*! \brief A location-aware function to add data to the beginning of a linked list.

    \param linkedList A linked list.

    \param data The data that we want to add to the beginning of linkedList.

    \return The link where data was placed in the linked list.
*/

Link* addFirst(LinkedList* linkedList, void* data)
{
    assert(linkedList != NULL);

    return addAfter(linkedList->head, data);
}

/*! \brief A location-aware function to add data to the end of a linked list.

    \param linkedList A linked list.

    \param data The data that we want to add to the end of linkedList.

    \return The link where data was placed in the linked list.
*/

Link* addLast(LinkedList* linkedList, void* data)
{
    assert(linkedList != NULL);

    return addBefore(linkedList->tail, data);
}

/*! \brief return the first piece of data in a linked list

    \param linkedList A linked list.

    \return The data in the first link of the linked list.
*/

void* getFirst(LinkedList* linkedList)
{
#ifdef DEBUG
    printf("getFirst...\n");
#endif
    assert(linkedList != NULL);
    assert(!isEmpty(linkedList));

    return linkedList->head->next->data;
}

/*! \brief remove the first link from a linked list

    \param linkedList A linked list.

    \return The first link of the linked list.
*/

Link* removeFirst(LinkedList* linkedList)
{
    assert(linkedList != NULL);

    if(!isEmpty(linkedList))
        return removeLink(linkedList->head->next);

    return NULL;
}

/*! \brief Remove and return the last link from a linked list.

    \param linkedList A linked list.

    \return The last link of the linked list.
*/

Link* removeLast(LinkedList* linkedList)
{
    assert(linkedList != NULL);

    if(!isEmpty(linkedList))
        return removeLink(linkedList->tail->prev);

    return NULL;
}

/*! \brief delete the last link in the linked list

    \param linkedList A linked list.
*/

void deleteLast(LinkedList* linkedList)
{
    assert(linkedList != NULL);
    if(!isEmpty(linkedList))
        delete(linkedList->tail->prev);

    return;
}

/*! \brief Print the first 10 items in the linked list

    \param linkedList A linked list.

    \param printFunc A function to print the data elements in
                     the linked list.
*/

void printListAbbv(LinkedList* linkedList, void (*printFunc)(void*))
{
#ifdef DEBUG
    printf("printListAbbv...\n");
#endif
    Link* curr = linkedList->head;
    curr = curr->next;

    int count = 0;

    while(!isTail(curr) && count != 10)
    {
        count++;
        printFunc(curr->data);
        if(!isTail(curr->next))
        {
            printf(",");
        }
        curr = curr->next;
    }

    if(!isTail(curr))
    {
        printf("... plus %d more", length(linkedList)-10);
    }

    printf("\n");
   
}

/*! \brief Print the items in the linked list.

    \param linkedList A linked list.

    \param printFunc A function to print the data elements in
                     the linked list.
*/

void printList(LinkedList* linkedList, void (*printFunc)(void*))
{
#ifdef DEBUG
    printf("printList...\n");
#endif
    Link* curr = linkedList->head->next;

    while(!isTail(curr))
    {
        printFunc(curr->data);
        if(!isTail(curr->next))
        {
            printf(" ");
        }
        curr = curr->next;
    }

    printf("\n");
   
}

/*! \brief Compute the number of data elements in a linked list.

    \param linkedList A linked list.

    \return The number of data elements in the linked list.
*/

int length(LinkedList* linkedList)
{
#ifdef DEBUG
    printf("length...\n");
#endif

#ifndef SMARTLENGTH
    int length = 0;
    Link* curr = linkedList->head->next;

    while(!isTail(curr))
    {
        length++;
        curr = curr->next;
    }

    return length;
#else
    assert(linkedList != NULL);
    return linkedList->length;
#endif
}

/*! \brief Determine if a linked list is empty.

    \param linkedList A linked list.

    \return Non-zero if the linked list is empty, zero otherwise.
*/

int isEmpty(LinkedList* linkedList)
{
#ifdef DEBUG
    printf("isEmpty...\n");
#endif

    assert(linkedList != NULL);

    return isTail(linkedList->head->next);
}

//@}
