#ifndef TASK_HPP
#define TASK_HPP

#include <vector>
#include <iostream>
#include <sstream>
#include <string>

#include "MPIPool.hpp"

using std::string;
using std::vector;
using std::stringstream;

class Task: public IMPITask {
public:
  Task(int v) :
      val(v) {
    ;
  }
  explicit Task(string s) {
    deserialise(s);
  }

  string run() {
    // just do some junk processing
    int junk;
    for (size_t i = 0; i < 1000000000; i++)
      junk += i;

    stringstream ss;
    ss << "foobar is " << this->val;
    return ss.str();
  }
  void deserialise(string msg) {
    stringstream ss;
    ss << msg;
    ss >> this->val;
  }
  string serialise() const {
    stringstream ss;
    ss << this->val;
    return ss.str();
  }
private:
  int val;
};

#endif
