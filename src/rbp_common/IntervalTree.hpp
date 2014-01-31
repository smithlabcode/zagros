/**
 \file IntervalTree.hpp
 \brief  An Interval Tree is a data structure for quickly determining
         the set of intervals that intersect a given point or
         interval. These classes are templates with two parameters, T
         and R. T is the type of the objects being stored in the tree
         and R is the type of their indices (i.e. int, or long would
         be a common choice). In theory the only restrictions are that
         R is ordinal (specifically, must define +, - and /). T can be
         any type as long as the user can provide pointers to functions
         'R getStart(T)' and 'R getEnd(T)' (trivial in most cases, and
         should be do-able with a functor.

 \authors Philip J. Uren

 \section copyright Copyright Details

 Copyright (C) 2011
 University of Southern California,
 Philip J. Uren

 \section license License Details

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 \section bugs Known Bugs

 \section history Revision History
   PJU -- 09 Mar 2011 -- Original port from Python implementation
   PJU -- 10 Mar 2011 -- Added template parameter R
   PJU -- 28 Jan 2014 -- Fixed some pass-by-value bugs that were slowing
                         this down; minor documentation clean up; other misc
                         optimizations to improve speed.
**/

#ifndef INTERVAL_TREE_HPP
#define INTERVAL_TREE_HPP

#include <vector>
#include <string>
#include <exception>
#include <sstream> 
#include <iostream>
#include <assert.h>

/******************************************************************************
 * Class definitions and prototypes
 *****************************************************************************/

/**
 * \summary Exception class for interval trees
 */
class IntervalTreeError: public std::exception
{
public:
	IntervalTreeError (const char* msg){this->msg = msg;}
  virtual const char* what() const throw(){return this->msg;}
private:
  const char* msg;
};

/**
 * \summary Functor for use in sorting Intervals of type <T> using a
 *          comparison function provided at construction time
 */
template <class T, class R>
class IntervalComparator
{
public:
  IntervalComparator(R (*compFunc)(T)) { this->compFunc = compFunc; }
  bool operator()(const T &i1, const T &i2){
    return this->compFunc(i1) < this->compFunc(i2);
  }
private:
	R (*compFunc)(T);
};

/**
 * \summary This class represents a single vertex (node) within an interval
 *          tree, it stores a set of intervals sorted by start and end.
 */
template <class T, class R>
class IntervalTreeNode {
public :
  IntervalTreeNode(const std::vector<T> &intervals, const double mid,
                   R (*getStart)(T), R (*getEnd)(T));
  IntervalTreeNode(const IntervalTreeNode &t);
  std::string toString();
  std::vector<T> starts;
  std::vector<T> ends;
  double mid;
private:
  R (*getStart)(T); 
  R (*getEnd)(T);
};

/**
 * \summary The actual IntervalTree class
 */
template <class T, class R>
class IntervalTree {
public:
  // constructors, destructors
  IntervalTree();
  IntervalTree(const std::vector<T> &intervals,
               R (*getStart)(T), R (*getEnd)(T));
  IntervalTree(const IntervalTree &r);
  ~IntervalTree();
  
  // inspectors
  void intersectingPoint(const R &point, std::vector<T> &v) const;
  void intersectingInterval(const R &start, const R &end, std::vector<T> &v) const;
  void squash(std::vector<T> &v) const;
  const size_t size() const;
  const std::string toString() const;
private:
  size_t treeSize;
  IntervalTreeNode<T,R>* data;
  IntervalTree<T,R>* left;
  IntervalTree<T,R>* right; 
  R (*getStart)(T); 
  R (*getEnd)(T);
};

/******************************************************************************
 * IntervalTreeNode class implementation
 *****************************************************************************/

/****
 * \summary IntervalTreeNode constructor
 * \param intervals the intervals to store at this node. Copies of the T's in
 *                  this vector are made, and these are what the tree stores.
 * \param mid       TODO
 * \param getStart  TODO
 * \param getEnd    TODO
 * TODO some of this can be done as initializations, rather than in the body
 *      of the c'tor.
 */
template <class T, class R> 
IntervalTreeNode<T,R>::IntervalTreeNode(const std::vector<T> &intervals,
                                        const double mid, R (*getStart)(T),
                                        R (*getEnd)(T)) {
  this->getStart = getStart;
  this->getEnd = getEnd;
  IntervalComparator<T,R> startComp = IntervalComparator<T,R>(getStart);
  IntervalComparator<T,R> endComp = IntervalComparator<T,R>(getEnd);
  this->starts = std::vector<T> (intervals);
  sort (this->starts.begin(), this->starts.end(), startComp);
  this->ends = std::vector<T> (intervals);
  sort (this->ends.begin(), this->ends.end(), endComp);
  this->mid = mid;
}

/**
 * \brief IntervalTreeNode copy constructor
 */
template <class T, class R>
IntervalTreeNode<T,R>::IntervalTreeNode(const IntervalTreeNode &r) {
  this->starts = std::vector<T>(r.starts);
  this->ends = std::vector<T>(r.ends);
  this->mid = r.mid;
  this->getStart = r.getStart;
  this->getEnd = r.getEnd;
}

/**
 * \summary construct a string representation of an IntervalTreeNode
 * \return the string representation
 */
template <class T, class R> 
std::string 
IntervalTreeNode<T,R>::toString() {
  std::ostringstream s;
  s << "mid: " << this->mid;
  s << std::endl;
  s << "intervals sorted by start:" << std::endl;
  for (size_t i = 0; i < this->starts.size(); ++i) {
    s << "(" << this->getStart(starts[i]) << " - " << this->getEnd(starts[i])
      << ")";
    s << std::endl;
  }
  s << "intervals sorted by end:" << std::endl;
  for (size_t i = 0; i < this->ends.size(); ++i) {
    s << "(" << this->getStart(ends[i]) << " - " << this->getEnd(ends[i])
      << ")";
    s << std::endl;
  }
  return s.str();
}

/******************************************************************************
 * IntervalTree class implementation
 *****************************************************************************/


/**
 * \summary: Constructor for IntervalTree.
 * \param intervals list of intervals, doesn't need to be sorted in any way.
 * \param getStart  TODO
 * \param getEnd    TODO
 * \throws IntervalTreeError: if no intervals are provided (NULL or empty)
 * \todo move some things into initializer list
 */
template <class T, class R> 
IntervalTree<T,R>::IntervalTree(const std::vector<T> &intervals,
								R (*getStart)(T), R (*getEnd)(T)) :
		treeSize(intervals.size()), left(NULL), right(NULL), getStart(getStart),
		getEnd(getEnd) {
  // can't build a tree with no intervals...
  if (intervals.size() <= 0)
    throw IntervalTreeError("Interval tree constructor got empty set of " 
														"intervals");
  
  // we need find the median; it's not easily done without re-ordering the
  // list, but we want to keep <intervals> const, so make a copy.
  // TODO -- PJU: there must be a way to do this efficiently without copying
  size_t midIndex = intervals.size() / 2;
  std::vector<T> cIntervals = std::vector<T>(intervals);
  IntervalComparator<T,R> startComp = IntervalComparator<T,R>(getStart);
  std::nth_element (cIntervals.begin(), cIntervals.begin()+midIndex,
                    cIntervals.end(), startComp);
  T midInt = cIntervals[midIndex];
  R mid = ((this->getEnd(midInt) - this->getStart(midInt)) / 2) 
  						 + this->getStart(midInt);
  
  std::vector<T> here, lt, rt;
  for (size_t i = 0; i < cIntervals.size(); ++i) {
    // place all intervals that end before <mid> into the left subtree
    if (this->getEnd(cIntervals[i]) < mid) lt.push_back(cIntervals[i]);
    else {
      // place all intervals that begin after <mid> into the right subtree
      if (this->getStart(cIntervals[i]) > mid) rt.push_back(cIntervals[i]);
        // everything else must overlap mid, so we keep it here
        else here.push_back(cIntervals[i]);
    }
  }

  // sanity check -- the mid point must at least overlap the middle interval
  if (here.size() <= 0) {
  	std::ostringstream msg;
  	msg << "fatal error: picked mid point at " << mid << " but this failed "
  	    << "to intersect anything!";
    throw IntervalTreeError(msg.str().c_str());
  }
      
  if (lt.size() > 0) this->left = new IntervalTree(lt, getStart, getEnd);
  if (rt.size() > 0) this->right = new IntervalTree(rt, getStart, getEnd);
  this->data = new IntervalTreeNode<T,R>(here, mid, this->getStart,
                                         this->getEnd);
}

/**
 * \brief IntervalTree copy constructor
 */
template <class T, class R>
IntervalTree<T,R>::IntervalTree(const IntervalTree &r) : data(NULL),
                                                         left(NULL),
                                                         right(NULL) {
  this->treeSize = r.treeSize;
  if (r.data != NULL) this->data = new IntervalTreeNode<T,R>(*(r.data));
  else this->data = NULL;
  if (r.left != NULL) this->left = new IntervalTree<T,R>(*(r.left));
  if (r.right != NULL) this->right = new IntervalTree<T,R>(*(r.right));
  this->getStart = r.getStart;
  this->getEnd = r.getEnd;
}

/**
 * \summary Destructor for IntervalTree
 */
template <class T, class R> 
IntervalTree<T,R>::~IntervalTree() {
	delete this->data;
	if (this->left != NULL) delete this->left;
	if (this->right != NULL) delete this->right;
}

/**
 * \summary given a point, determine which set of intervals in the tree
 *          are intersected.
 * \param point         the point of intersection to test against
 * \param intersect     intersecting intervals will be added to this vector;
 *                      any existing elements will NOT be cleared.
 */
template <class T, class R> 
void
IntervalTree<T,R>::intersectingPoint(const R &point,
                                     std::vector<T> &intersect) const {
	// perfect match, everything here overlaps
	if (point == this->data->mid) {
	  intersect.insert(this->data->ends.begin(), this->data->ends.end());
	}
	else if (point > this->data->mid) {
      // we know all intervals in this->data begin before p (if they began
      // after p, they would have not included mid) we just need to find
      // those that end after p
      for (size_t i = 0; i < this->data->ends.size(); ++i) {
        if (this->getEnd(this->data->ends[i]) >= point)
          intersect.push_back(this->data->ends[i]);
        else
          break;
      }
      // recurse down the right (after mid) branch
      if (this->right != NULL) this->right->intersectingPoint(point, intersect);
	} else { //point < this->data->mid
      // we know all intervals in this->data end after p (if they ended before
      // p, they would have not included mid) we just need to find those that
      // start before p
      for (size_t i = 0; i < this->data->starts.size(); ++i) {
        if (this->getStart(this->data->starts[i]) <= point)
          intersect.push_back(this->data->starts[i]);
        else
          break;
      }
      // recurse down the left (before mid) branch
      if (this->left != NULL) this->left->intersectingPoint(point, intersect);
	}
}

/**
 * \summary given an interval, determine which set of intervals in the tree
 *          are intersected. The result is a vector of const pointers to the
 *          intersected intervals; this provides fast access without needing
 *          to copy anything, but we can't allow the caller to modify them,
 *          as it may invalidate the tree.
 * \param start TODO
 * \param end   TODO
 * \param intersect     intersecting intervals will be added to this vector;
 *                      any existing elements will NOT be cleared.
 */
template <class T, class R> 
void
IntervalTree<T,R>::intersectingInterval(const R &start,	const R &end,
                                        std::vector<T> &intersect) const {
  // find all intervals in this node that intersect start and end
  for (size_t i = 0; i < this->data->starts.size(); ++i) {
    T &current = this->data->starts[i];
    if ((this->getEnd(current) > start) && (this->getStart(current) < end)) {
      intersect.push_back(current);
    }
  }
	
  // process left subtree (if we have one) if the requested 
  // interval begins before mid
  if ((this->left != NULL) && (start <= this->data->mid))
  	this->left->intersectingInterval(start, end, intersect);

  // process right subtree (if we have one) if the requested interval 
  // ends after mid
  if ((this->right != NULL) && (end >= this->data->mid))
  	this->right->intersectingInterval(start, end, intersect);
}

/**
 * \summary squash the tree; i.e. populate a vector with all of the intervals
 *          in the tree. This is not destructive, the original tree remains and
 *          the resultant vector contains copies of the elements in the tree.
 *          This might be very expensive if copying a T is expensive.
 * \param   intervals   will contain all of the intervals from the tree after
 *                      the call; any data in the vector beforehand will NOT
 *                      be cleared
 */
template <class T, class R>
void
IntervalTree<T,R>::squash(std::vector<T> &intervals) const {
  intervals.insert(intervals.end(), this->data->starts.begin(),
                   this->data->starts.end());
  if (this->left != NULL) this->left->squash(intervals);
  if (this->right != NULL) this->right->squash(intervals);
}

/**
 * \summary get the size of the tree, defined as the number of intervals
 *          stored in the tree.
 * \return  the size of (number of intervals in) the tree
 */
template <class T, class R>
const size_t
IntervalTree<T,R>::size() const {
  return this->treeSize;
}

/**
 * \summary return a string representation of an IntervalTree
 */
template <class T, class R> 
const std::string 
IntervalTree<T,R>::toString() const {
  std::string left = "<EMPTY>", right = "<EMPTY>";
  if (this->left != NULL) left = this->left->toString();
  if (this->right != NULL) right = this->right->toString();
  return this->data->toString() + "\n** left ** " + left +
         "\n** right ** " + right;
}
	
#endif
