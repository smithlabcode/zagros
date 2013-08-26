// TODO what is this file for? Can we delete it?

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <float.h> 
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include  "part_func.hpp"

using std::string;
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;
using std::endl;
using std::vector;

int main(int argc, char *argv[]) {
  /*float energy = -1;

   if (argc > 1) {
   string s;
   s = (string)argv[1];
   cout << s << endl;
   int TempNumOne = s.length();
   char seq1[s.length()];
   for (int i = 0; i <= TempNumOne; i++) {
   seq1[i] = s[i];
   if (seq1[i] == 'N')
   return 0;
   }

   char *str1 = (char* ) malloc(sizeof(char)*(s.length()+1));

   s = (string)argv[2];
   cout << s << endl;
   TempNumOne = s.length();
   for (int i = 0; i <= TempNumOne; i++)
   str1[i] = s[i];
   float e1, kT;

   e1 = 0;
   MC mc;
   kT = ( mc.getTemperature() + 273.15)*1.98717/1000.0;  // kT in kcal/mol
   mc.init_pf_fold(s.length(), energy);
   e1 = mc.pf_fold(seq1, str1);
   cout << str1 << endl;
   mc.printMatrix("whole_big_probs.out", 2, s.length(), 0);
   vector<double> T;
   T.resize(s.length());
   mc.getProbVector(T);
   for ( int i = 0; i < s.length(); i++)
   cout << T[i] << ",";
   cout << endl;
   }*/
  return 0;
}
