#ifndef _DJS_LINKEDLIST_H_
#define _DJS_LINKEDLIST_H_

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

/*! \file LinkedList.h

    \brief see LinkedList.c

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

struct Link;

/*! \struct LinkedList

    \brief Stores a linked list, with sentinel links for head and tail.
           These sentinels contain dummy values. 

*/

struct LinkedList
{
    struct Link* head; //!< head of the linked list, a dummy sentinel
    struct Link* tail; //!< tail of the linked list, a dummy sentinel
#ifdef SMARTLENGTH
    int length; //!< the number of data items in the linked list, if we are maintaining it.
#endif
};

/*! \struct Link

    \brief Stores data and pointers to next and previous links.

*/

struct Link
{
    void* data; //!< arbitrary data stored in the link
    struct Link* next; //!< the previous link in the chain
    struct Link* prev; //!< the next link in the chain
#ifdef SMARTLENGTH
    struct LinkedList* linkedList; //!< the linked list that this link belongs to.
#endif
};

typedef struct Link Link;

typedef struct LinkedList LinkedList;

int isHead(Link* list);

int isTail(Link* list);

void* delete(Link* list);

Link* addAfter(Link* list, void* data);

Link* addBefore(Link* list, void* data);

void addLinkBefore(Link* list, Link* newLink);

Link* removeLink(Link* list);

LinkedList* createLinkedList(void);

void destroyLinkedList(LinkedList* linkedList);

void copyLinkedList(LinkedList* destination, 
                    LinkedList* source);

int contains(LinkedList* linkedList, void* data, int (*comparator)(void*,void*));

int equal( LinkedList* list1, 
           LinkedList* list2, 
           int (*comparator)(void*,void*));

void restoreLinksWithReferences(LinkedList* list);

Link* addFirst(LinkedList* linkedList, void* data);

Link* addLast(LinkedList* linkedList, void* data);

Link* removeFirst(LinkedList* linkedList);

Link* removeLast(LinkedList* linkedList);

void deleteLast(LinkedList* linkedList);

void* getFirst(LinkedList* linkedList);

void printListAbbv(LinkedList* linkedList, void (*printFunc)(void*));

void printList(LinkedList* linkedList, void (*printFunc)(void*));

int length(LinkedList* linkedList);

int isEmpty(LinkedList* linkedList);

#endif
