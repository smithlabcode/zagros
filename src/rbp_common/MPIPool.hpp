/**
 \file MPIPool.hpp

 \brief This header declares a mechanism to coordinate multiple processes into
 a compute pool using MPI.

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

#ifndef MPIPOOL_H
#define MPIPOOL_H

#include <vector>
#include <queue>
#include <string>
#include "mpi.h"

/**
 * \brief An interface that defines what public members a class must have
 *        to be submitted to an MPIPool for parallel distributed execution
 */
class IMPITask {
public:
  /** \brief derived must implement destructor to stop destruction via base **/
  virtual ~IMPITask() {
    ;
  }
  /** \brief execute this task and return the result as a string **/
  virtual std::string run() = 0;
  /** \brief take a serialised string of this class and populate the object **/
  virtual void deserialise(std::string msg) = 0;
  /** \brief produce a serialised representation of this object **/
  virtual std::string serialise() const = 0;
};

/**
 * \brief A class that coordinates communication between processes submitting
 *        jobs and those processing them
 *
 * This class has several important restrictions that stem from the fact that
 * it was written to solve a very specific problem quickly and without undue
 * complexity that, while it would be nice in the general case, was not needed
 * in the specific one. Firstly, it requires all MPI processes in the current
 * MPI program to join the pool as workers except for the rank 0 process; if
 * any fail to do so, or unexpectedly die, the process submitting jobs will
 * wait indefinitely for them. Secondly, the only process permitted to submit
 * jobs to the pool is the rank 0 process. Finally, and importantly, this class
 * is <strong>not</strong> thread safe.
 */
class MPIPoolQueue {
public:
  /*** Constructors, destructors and object initialisation ***/
  MPIPoolQueue();

  /*** Defining the pool interface ***/
  /*** Inspectors ***/
  bool isIdle() const;
  bool isDone() const;
  /*** Mutators ***/
  void join();
  void despatchJob(IMPITask const* const t);
  void wait();
  void release();

  /*** defining the the queue-like interface for job results ***/
  bool isEmpty() const;
  std::string dequeue();

private:
  /*** Private Types ***/
  /** \brief Tags for separating messages into control and data **/
  enum MessageTag {
    CONTROL, DATA
  };
  /** \brief Predefined control messages **/
  enum ControlMessages {
    DATA_MSG_PENDING, FORCE_QUIT, NOTIFY_QUIT
  };

  /*** Private member functions ***/
  void listen() const;
  void receive(int rank);
  int makeCapacity();

  /*** Private member variables ***/
  /** \brief True if the pool is ready to receive jobs, false otherwise **/
  bool poolAvailable;
  /** \brief MPI_Request objects for storing the result of posted IRecvs **/
  std::vector<MPI_Request> requests;
  /** \brief Buffers for control messages indicating job completion **/
  std::vector<int> requestCompletionMessages;
  /** \brief Booleans indicating which worker processes are busy **/
  std::vector<bool> freeProcesses;
  /** \brief Queue for storing the results received from worker processes **/
  std::queue<std::string> outputQueue;
};

#endif
