/**
 \file MPIPool.cpp

 \brief This source file defines the details of a mechanism to coordinate
 multiple processes into a compute pool using MPI.

 \author Philip J. Uren

 \section copyright Copyright Details
 Copyright (C) 2012
 University of Southern California,
 Philip J. Uren, Emad Bahrami Samani, Andrew D. Smith

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

 \section revisions Revision History
 */

#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mpi.h"

#include "MPIPool.hpp"
#include "Task.hpp"

using std::string;
using std::vector;
using std::stringstream;
using std::cerr;
using std::endl;

/**
 * \brief Construct an MPIPoolQueue that accepts jobs from rank 0 and
 *        distributes them to all other ranks.
 */
MPIPoolQueue::MPIPoolQueue() :
    poolAvailable(true) {
  int p;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  requests = vector<MPI_Request>(p, MPI_Request());
  requestCompletionMessages = vector<int>(p, 0); // TODO magic number
  freeProcesses = vector<bool>(p, true);
  freeProcesses[0] = false; // this is rank 0, it's never free
}

/**
 * \brief Returns true if the pool has no pending work to do
 */
bool MPIPoolQueue::isIdle() const {
  for (size_t j = 1; j < freeProcesses.size(); j++) {
    if (!freeProcesses[j])
      return false;
  }
  return true;
}

/**
 * \brief Returns true if the pool is idle and its output queue is empty
 */
bool MPIPoolQueue::isDone() const {
  return (isIdle() && isEmpty());
}

/**
 * \brief Returns true if there are no completed jobs waiting in the
 *        pool's output queue
 */
bool MPIPoolQueue::isEmpty() const {
  return this->outputQueue.size() == 0;
}

/**
 * \brief Remove and return the first finished job from the queue
 * \todo  throw exception if the queue is empty
 */
string MPIPoolQueue::dequeue() {
  string frnt = this->outputQueue.front();
  this->outputQueue.pop();
  return frnt;
}

/**
 * \brief  Wait for at least one worker to become free
 * \return Rank of free worker process
 *
 * After a call to this function, the pool is guaranteed to have at least one
 * free process available for work, the rank of which is returned as an int
 */
int MPIPoolQueue::makeCapacity() {
  for (size_t j = 1; j < freeProcesses.size(); j++) {
    if (freeProcesses[j])
      return j;
  }

  // Nothing is free, block until something finishes. We have to make this
  // temporary array, passing the address of the first item in the vector
  // fails
  MPI_Request tmp[this->requests.size() - 1];
  for (size_t i = 1; i < this->freeProcesses.size(); i++)
    tmp[i - 1] = this->requests[i];
  int finIndx = -1;
  MPI_Waitany(this->freeProcesses.size() - 1, tmp, &finIndx, MPI_STATUS_IGNORE);
  finIndx += 1; // adjust for rank 0

  // Now a worker has signaled that it's finished it's work...
  this->receive(finIndx);
  return finIndx;
}

/**
 * \brief Wait for all outstanding jobs to complete and their results to be
 *        placed into the output queue.
 * \note  This is a blocking operation
 */
void MPIPoolQueue::wait() {
  for (size_t j = 1; j < freeProcesses.size(); j++) {
    if (!freeProcesses[j]) {
      MPI_Wait(&(requests[j]), MPI_STATUS_IGNORE);
      this->receive(j);
    }
  }
}

/**
 * \brief   Receive the data message from process <rank>
 * \note:   This will block until the data message is ready. It does not
 *          receive the control message - assume that's already been received
 *
 * \todo    throw exception if rank is marked as being currently free for work
 */
void MPIPoolQueue::receive(int rank) {
  // probe for the size of the data message
  MPI_Status status;
  MPI_Probe(rank, DATA, MPI_COMM_WORLD, &status);
  int msgSize;
  MPI_Get_count(&status, MPI_CHAR, &msgSize);

  // make space for the data message and collect it
  char* buf = (char*) malloc(sizeof(char) * msgSize);
  MPI_Recv(
      buf, msgSize, MPI_CHAR, rank, DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  this->outputQueue.push(buf);
  this->freeProcesses[rank] = true;
}

/**
 * \brief Despatch the given task to a worker process for execution.
 *
 * This is a blocking operation, it won't return until the job has been
 * despatched to a free process. That is to say, there is no queuing of
 * incoming jobs. We don't keep any references or pointers to the given
 * concrete IMPITask object, so the caller is free to do whatever they want
 * with it afterwards
 */
void MPIPoolQueue::despatchJob(IMPITask const* const t) {
  if (!this->poolAvailable) {
    cerr << "oops, pool has already been closed, "
        << "can't despatch any more jobs, sorry" << endl;
    // TODO throw exception here
    return;
  }

  int freeServer = this->makeCapacity();
  int f = DATA_MSG_PENDING;
  MPI_Send(&f, 1, MPI_INT, freeServer, CONTROL, MPI_COMM_WORLD);
  string jobstr = t->serialise();
  char* jobStrC = new char[jobstr.size() + 1];
  strcpy(jobStrC, jobstr.c_str());
  MPI_Send(
      jobStrC, strlen(jobstr.c_str()) + 1, MPI_CHAR, freeServer, DATA,
      MPI_COMM_WORLD);
  delete jobStrC;

  // post that we're waiting for a message back and mark this server as busy
  int* msgResBuff = &(this->requestCompletionMessages[freeServer]);
  MPI_Request* msgRequest = &(this->requests[freeServer]);
  MPI_Irecv(
      msgResBuff, 1, MPI_INT, freeServer, CONTROL, MPI_COMM_WORLD, msgRequest);
  this->freeProcesses[freeServer] = false;
  cerr << "despatched job to " << freeServer << endl;
}

/**
 * \brief Join this pool as a worker process
 *
 * The calling process joins the pool as a worker process and will not return
 * from this call until until it receives a control message indicating it
 * should quit the pool.
 */
void MPIPoolQueue::join() {
  this->listen();
}

/**
 * \brief Release all worker processes from the pool. The pool cannot be used
 *        used any further after doing this.
 */
void MPIPoolQueue::release() {
  for (size_t i = 1; i < this->freeProcesses.size(); i++) {
    int f = FORCE_QUIT;
    MPI_Send(&f, 1, MPI_INT, int(i), CONTROL, MPI_COMM_WORLD);
  }
  this->poolAvailable = false;
}

/**
 * \brief Listen for incoming jobs, process them and return the results
 *
 *        The calling process is put to work accepting incoming jobs,
 *        performing them and returning the results to the process that
 *        submitted the job.
 */
void MPIPoolQueue::listen() const {
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  while (true) {
    // wait for control code from despatcher -- for now we just always know
    // the despatcher is at rank 0. This call blocks.
    int controlMsg;
    MPI_Recv(
        &controlMsg, 1, MPI_INT, 0, CONTROL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (controlMsg == FORCE_QUIT)
      break;
    else if (controlMsg == DATA_MSG_PENDING) {
      // probe for the message from the client, again for now we just know
      // the client is at rank 0
      MPI_Status status;
      MPI_Probe(0, DATA, MPI_COMM_WORLD, &status);
      // how big is the message?
      int msgSize;
      MPI_Get_count(&status, MPI_CHAR, &msgSize);
      // Allocate a buffer and collect the message
      char* buf = new char[msgSize + 1];
      MPI_Recv(
          buf, msgSize, MPI_CHAR, 0, DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // do the job
      // TODO add another message which gives type; allow types to be
      // registered
      Task tsk = Task(string(buf));
      string res = tsk.run();

      // send back the result -- we have to fudge f because MPI_Send wants the
      // message as void pointer, who knows why...
      int f = DATA_MSG_PENDING;
      MPI_Send(&f, 1, MPI_INT, 0, CONTROL, MPI_COMM_WORLD);
      char* cres = new char[res.size() + 1];
      strcpy(cres, res.c_str());
      MPI_Send(
          cres, strlen(res.c_str()) + 1, MPI_CHAR, 0, DATA, MPI_COMM_WORLD);

      // cleanup after ourselves
      delete cres;
      delete buf;
    }
  }
}
