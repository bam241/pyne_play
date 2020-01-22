// This file is composed of the following original files:

//   /home/mouginot/work/app/pyne/license.txt
//   /home/mouginot/work/app/pyne/src/utils.cpp
//   /home/mouginot/work/app/pyne/src/state_map.cpp
//   /home/mouginot/work/app/pyne/src/nucname.cpp
//   /home/mouginot/work/app/pyne/src/jsoncpp.cpp
//   /home/mouginot/work/app/pyne/src/jsoncustomwriter.cpp
//   /home/mouginot/work/app/pyne/src/material.cpp
//   /home/mouginot/work/app/pyne/src/material_library.cpp
//   /home/mouginot/work/app/pyne/src/particle.cpp
//   /home/mouginot/work/app/pyne/src/tally.cpp

// PyNE amalgated source http://pyne.io/
#include "pyne.h"

//
// start of /home/mouginot/work/app/pyne/license.txt
//
// Copyright 2011-2019, the PyNE Development Team. All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
// 
//    1. Redistributions of source code must retain the above copyright notice, this list of
//       conditions and the following disclaimer.
// 
//    2. Redistributions in binary form must reproduce the above copyright notice, this list
//       of conditions and the following disclaimer in the documentation and/or other materials
//       provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE PYNE DEVELOPMENT TEAM ``AS IS'' AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// The views and conclusions contained in the software and documentation are those of the
// authors and should not be interpreted as representing official policies, either expressed
// or implied, of the stakeholders of the PyNE project or the employers of PyNE developers.
// 
// -------------------------------------------------------------------------------
// The files cpp/measure.cpp and cpp/measure.hpp are covered by:
// 
// Copyright 2004 Sandia Corporation.  Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
// retains certain rights in this software.
// 
// https://press3.mcs.anl.gov/sigma/moab-library
// 
// -------------------------------------------------------------------------------
// The files in fortranformat/ are covered by: 
// 
// The MIT License. Copyright (c) 2011 Brendan Arnold
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// 
// https://bitbucket.org/brendanarnold/py-fortranformat/src/
// //
// end of /home/mouginot/work/app/pyne/license.txt
//


//
// start of /home/mouginot/work/app/pyne/src/utils.cpp
//
// General Library
#ifndef PYNE_IS_AMALGAMATED
extern "C" double endftod_(char *str, int len);
#endif
#include <iomanip>

#ifndef PYNE_IS_AMALGAMATED
#include "utils.h"
#endif


// PyNE Globals

std::string pyne::PYNE_DATA = "";
std::string pyne::NUC_DATA_PATH = "";
std::string pyne::VERSION = "0.5.11";

void pyne::pyne_start() {
#if defined __WIN_MSVC__
  char * tmpPYNE_DATA;
  size_t lenPYNE_DATA;
  errno_t errPYNE_DATA = _dupenv_s(&tmpPYNE_DATA, &lenPYNE_DATA, "PYNE_DATA");
  if (errPYNE_DATA)
    tmpPYNE_DATA = (char *) "<NOT_FOUND>";
  PYNE_DATA = (std::string) tmpPYNE_DATA;

  char * tmpNUC_DATA_PATH;
  size_t lenNUC_DATA_PATH;
  errno_t errNUC_DATA_PATH = _dupenv_s(&tmpNUC_DATA_PATH, &lenNUC_DATA_PATH, "NUC_DATA_PATH");
  if (errPYNE_DATA)
    tmpNUC_DATA_PATH = (char *) "<NOT_FOUND>";
  NUC_DATA_PATH = (std::string) tmpNUC_DATA_PATH;
#else
  char * tmppath;
  tmppath = getenv("PYNE_DATA");
  if (tmppath == NULL)
      tmppath = (char *) "<NOT_FOUND>";
  PYNE_DATA = std::string(tmppath);

  tmppath = getenv("NUC_DATA_PATH");
  if (tmppath == NULL)
      tmppath = (char *) "<NOT_FOUND>";
  NUC_DATA_PATH = std::string(tmppath);
#endif
  return;
}



// String Transformations
std::string pyne::to_str(int t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(unsigned int t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(double t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string pyne::to_str(bool t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}


int pyne::to_int(std::string s) {
  return atoi( s.c_str() );
}

double pyne::to_dbl(std::string s) {
  return strtod( s.c_str(), NULL );
}

double pyne::endftod_cpp(char * s) {
  // Converts string from ENDF only handles "E-less" format but is 5x faster
  int pos, mant, exp;
  double v, dbl_exp;

  mant = exp = 0;
  if (s[2] == '.') {
    // Convert an ENDF float
    if (s[9] == '+' or s[9] == '-') {
      // All these factors of ten are from place values.
      mant = s[8] + 10 * s[7] + 100 * s[6] + 1000 * s[5] + 10000 * s[4] + \
             100000 * s[3] + 1000000 * s[1] - 1111111 * '0';
      exp = s[10] - '0';
      // Make the right power of 10.
      dbl_exp = exp & 01? 10.: 1;
      dbl_exp *= (exp >>= 1) & 01? 100.: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e4: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e8: 1;
      // Adjust for powers of ten from treating mantissa as an integer.
      dbl_exp = (s[9] == '-'? 1/dbl_exp: dbl_exp) * 1.0e-6;
      // Get mantissa sign, apply exponent.
      v = mant * (s[0] == '-'? -1: 1) * dbl_exp;
    }
    else {
      mant = s[7] + 10 * s[6] + 100 * s[5] + 1000 * s[4] + 10000 * s[3] + \
             100000 * s[1] - 111111 * '0';
      exp = s[10] + 10 * s[9] - 11 * '0';
      dbl_exp = exp & 01? 10.: 1;
      dbl_exp *= (exp >>= 1) & 01? 100.: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e4: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e8: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e16: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e32: 1;
      dbl_exp *= (exp >>= 1) & 01? 1.0e64: 1;
      dbl_exp = (s[8] == '-'? 1/dbl_exp: dbl_exp) * 1.0e-5;
      v = mant * (s[0] == '-'? -1: 1) * dbl_exp;
    }
  }

  // Convert an ENDF int to float; we start from the last char in the field and
  // move forward until we hit a non-digit.
  else {
    v = 0;
    mant = 1; // Here we use mant for the place value about to be read in.
    pos = 10;
    while (s[pos] != '-' and s[pos] != '+' and s[pos] != ' ' and pos > 0) {
      v += mant * (s[pos] - '0');
      mant *= 10;
      pos--;
    }
    v *= (s[pos] == '-'? -1: 1);
  }
  return v;
}

double pyne::endftod_f(char * s) {
#ifdef PYNE_IS_AMALGAMATED
  return endftod_cpp(s);
#else
  return endftod_(s, 12);
#endif
}

double (*pyne::endftod)(char * s) = &pyne::endftod_f;

void pyne::use_fast_endftod() {
  pyne::endftod = &pyne::endftod_cpp;
}

std::string pyne::to_upper(std::string s) {
  // change each element of the string to upper case.
  for(unsigned int i = 0; i < s.length(); i++)
    s[i] = toupper(s[i]);
  return s;
}

std::string pyne::to_lower(std::string s) {
  // change each element of the string to lower case
  for(unsigned int i = 0; i < s.length(); i++)
    s[i] = tolower(s[i]);
  return s;
}

std::ostringstream pyne::comment_line_wrapping(std::string line,
                                               std::string comment_prefix,
                                               int line_length) {
  std::ostringstream oss;

  line_length -= comment_prefix.length();
    
  // Include as is if short enough
  while (line.length() > line_length) {
    oss << comment_prefix << line.substr(0, line_length) << std::endl;
    line.erase(0, line_length);
  }

  if (line.length() > 0) {
    oss << comment_prefix << line << std::endl;
  }

  return oss;
}

std::string pyne::capitalize(std::string s) {
  unsigned int slen = s.length();
  if (slen == 0)
    return s;
  // uppercase the first character
  s[0] = toupper(s[0]);
  // change each subsequent element of the string to lower case
  for(unsigned int i = 1; i < slen; i++)
    s[i] = tolower(s[i]);
  return s;
}


std::string pyne::get_flag(char line[], int max_l) {
  char tempflag [10];
  for (int i = 0; i < max_l; i++)
  {
    if (line[i] == '\t' || line[i] == '\n' || line[i] == ' ' || line[i] == '\0')
    {
      tempflag[i] = '\0';
      break;
    }
    else
      tempflag[i] = line[i];
  }
  return std::string (tempflag);
}



std::string pyne::remove_substring(std::string s, std::string substr) {
  // Removes a substring from the string s
  int n_found = s.find(substr);
  while ( 0 <= n_found ) {
    s.erase( n_found , substr.length() );
    n_found = s.find(substr);
  }
  return s;
}


std::string pyne::remove_characters(std::string s, std::string chars) {
  // Removes all characters in the string chars from the string s
  for (int i = 0; i < chars.length(); i++ ) {
    s = remove_substring(s, chars.substr(i, 1) );
  }
  return s;
}


std::string pyne::replace_all_substrings(std::string s, std::string substr, std::string repstr) {
  // Replaces all instance of substr in s with the string repstr
  int n_found = s.find(substr);
  while ( 0 <= n_found ) {
    s.replace( n_found , substr.length(), repstr );
    n_found = s.find(substr);
  }
  return s;
}



std::string pyne::last_char(std::string s) {
    // Returns the last character in a string.
    return s.substr(s.length()-1, 1);
}


std::string pyne::slice_from_end(std::string s, int n, int l) {
  // Returns the slice of a string using negative indices.
  return s.substr(s.length()+n, l);
}


bool pyne::ternary_ge(int a, int b, int c) {
  // Returns true id a <= b <= c and flase otherwise.
  return (a <= b && b <= c);
}


bool pyne::contains_substring(std::string s, std::string substr) {
  // Returns a boolean based on if the sub is in s.
  int n = s.find(substr);
  return ( 0 <= n && n < s.length() );
}


std::string pyne::natural_naming(std::string name) {
  // Calculates a version on the string name that is a valid
  // variable name, ie it uses only word characters.
  std::string nat_name (name);

  // Replace Whitespace characters with underscores
  nat_name = pyne::replace_all_substrings(nat_name, " ",  "_");
  nat_name = pyne::replace_all_substrings(nat_name, "\t", "_");
  nat_name = pyne::replace_all_substrings(nat_name, "\n", "_");

  // Remove non-word characters
  int n = 0;
  while ( n < nat_name.length() ) {
    if ( pyne::words.find(nat_name[n]) == std::string::npos )
      nat_name.erase(n, 1);
    else
      n++;
  }

  // Make sure that the name in non-empty before continuing
  if (nat_name.length() == 0)
    return nat_name;

  // Make sure that the name doesn't begin with a number.
  if ( pyne::digits.find(nat_name[0]) != std::string::npos)
    nat_name.insert(0, "_");

  return nat_name;
}


std::vector<std::string> pyne::split_string(std::string particles_list, std::string delimiter) {
  std::vector<std::string> output_vector;
  size_t prev_pos = 0; //item start position
  size_t pos = 0; //item end position
 
  while( (pos = particles_list.find(delimiter, prev_pos)) != std::string::npos){
    output_vector.push_back(particles_list.substr(prev_pos, pos));
    prev_pos = pos + delimiter.length();
  }
  // catch list with a single particle
  if (pos == std::string::npos && prev_pos == 0 && particles_list.length() >0)
    output_vector.push_back(particles_list);

  return output_vector;
}



template<typename T>
std::string pyne::join_to_string(std::vector<T> vect, std::string delimiter){
  std::stringstream out;
  out << std::setiosflags(std::ios::fixed) << std::setprecision(6);
  
  // ensure there is at least 1 element in the vector
  if (vect.size() == 0)
    return out.str();
  // no delimiter needed before the first element
  out << vect[0];
  for( int i = 1; i < vect.size(); i++)
    out << delimiter << vect[i];
  return out.str();
}
template std::string pyne::join_to_string(std::vector<int> vect, std::string delimiter);
template std::string pyne::join_to_string(std::vector<double> vect,
                                 std::string delimiter);
template std::string pyne::join_to_string(std::vector<std::string> vect, std::string delimiter);

//
// Math Helpers
//

double pyne::slope(double x2, double y2, double x1, double y1) {
  // Finds the slope of a line.
  return (y2 - y1) / (x2 - x1);
}


double pyne::solve_line(double x, double x2, double y2, double x1, double y1) {
  return (slope(x2,y2,x1,y1) * (x - x2)) + y2;
}


double pyne::tanh(double x) {
  return std::tanh(x);
}

double pyne::coth(double x) {
  return 1.0 / std::tanh(x);
}



// File Helpers

bool pyne::file_exists(std::string strfilename) {
  // Thank you intarwebz for this function!
  // Sepcifically: http://www.techbytes.ca/techbyte103.html
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strfilename.c_str(), &stFileInfo);

  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  }
  else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }

  return(blnReturn);
}

// convert convert a filename into path+filename (for pyne)
std::string pyne::get_full_filepath(char* filename) {
  std::string file(filename);
  return pyne::get_full_filepath(file);
}

// convert convert a filename into path+filename (for pyne)
std::string pyne::get_full_filepath(std::string filename) {
  // remove all extra whitespace
  filename = pyne::remove_characters(" " , filename);
  // use stdlib call
  const char* full_filepath = realpath(filename.c_str(), NULL);
  return std::string(full_filepath);
}

// Message Helpers

bool pyne::USE_WARNINGS = true;

bool pyne::toggle_warnings(){
  USE_WARNINGS = !USE_WARNINGS;
  return USE_WARNINGS;
}

void pyne::warning(std::string s){
  // Prints a warning message
  if (USE_WARNINGS){
    std::cout << "\033[1;33m WARNING: \033[0m" << s << "\n";
  }
}




//
// end of /home/mouginot/work/app/pyne/src/utils.cpp
//


//
// start of /home/mouginot/work/app/pyne/src/state_map.cpp
//
//Mapping file for state ids to nuc ids
//This File was autogenerated!!
#ifndef PYNE_4HFU6PUEQJB3ZJ4UIFLVU4SPCM
#define PYNE_4HFU6PUEQJB3ZJ4UIFLVU4SPCM
namespace pyne {
namespace nucname {
#define TOTAL_STATE_MAPS 1016
std::map<int, int> state_id_map;
int map_nuc_ids [TOTAL_STATE_MAPS] = {110240001,
130240001,
90260001,
130260001,
130320002,
170340001,
170380001,
190380001,
190380015,
210420002,
210430001,
210440004,
230440001,
210450001,
210460002,
230460001,
210500001,
250500001,
250520001,
260520042,
260530022,
270530003,
270540001,
210560001,
210560004,
250580001,
270580001,
270580002,
230600000,
230600001,
250600001,
270600001,
260610004,
250620001,
270620001,
230640001,
250640002,
260650003,
260670002,
290670023,
290680003,
280690001,
280690008,
300690001,
340690004,
270700001,
290700001,
290700003,
350700006,
280710002,
300710001,
320710002,
310720002,
350720001,
300730001,
300730002,
320730002,
340730001,
360730004,
310740002,
350740002,
290750001,
290750002,
300750001,
320750002,
330750004,
280760004,
290760001,
350760002,
300770002,
320770001,
330770004,
340770001,
350770001,
300780004,
310780004,
350780004,
370780003,
390780001,
300790002,
320790001,
330790007,
340790001,
350790001,
360790001,
310800001,
350800002,
390800001,
390800003,
320810001,
340810001,
360810002,
370810001,
330820001,
340820015,
350820001,
370820001,
410820003,
340830001,
360830002,
370830002,
380830002,
390830001,
310840001,
350840001,
360840019,
360840061,
370840002,
390840002,
410840007,
360850001,
370850003,
380850002,
390850001,
400850002,
410850003,
410850005,
370860002,
380860014,
390860002,
410860001,
410860002,
380870001,
390870001,
400870002,
410870001,
430870005,
350880003,
410880001,
430880000,
430880001,
390890001,
400890001,
410890001,
420890002,
430890001,
370900001,
390900002,
400900003,
410900002,
410900007,
430900001,
430900006,
390910001,
400910040,
410910001,
420910001,
430910001,
440910001,
450910001,
410920001,
450920001,
390930002,
410930001,
420930016,
430930001,
440930001,
410940001,
430940001,
470940001,
470940002,
410950001,
430950001,
450950001,
460950005,
470950002,
390960005,
430960001,
450960001,
370970002,
390970001,
390970029,
410970001,
430970001,
450970001,
370980001,
390980005,
410980001,
450980001,
410990001,
430990002,
450990001,
470990002,
371000001,
391000004,
411000001,
411000009,
411000012,
431000002,
431000004,
451000004,
471000001,
451010001,
471010002,
411020001,
431020001,
451020005,
471020001,
441030005,
451030001,
471030002,
491030001,
411040004,
451040003,
471040001,
491040003,
451050001,
471050001,
491050001,
411060003,
451060001,
471060001,
491060001,
431070000,
461070002,
471070001,
491070001,
401080009,
411080003,
451080004,
471080002,
491080001,
461090002,
471090001,
491090001,
491090026,
451100001,
471100002,
491100001,
421110001,
461110002,
471110001,
481110003,
491110001,
491120001,
491120004,
491120009,
461130002,
471130001,
481130001,
491130001,
501130001,
451140005,
491140001,
491140005,
531140005,
461150001,
471150001,
481150001,
491150001,
521150001,
451160001,
471160001,
471160004,
491160001,
491160004,
511160003,
551160001,
441170003,
471170001,
481170002,
491170001,
501170002,
521170003,
471180004,
491180001,
491180003,
511180007,
531180002,
551180001,
441190002,
471190000,
471190001,
481190002,
491190001,
501190002,
511190072,
521190002,
551190001,
451200002,
471200002,
491200001,
491200002,
511200001,
531200013,
551200001,
571200000,
461210001,
481210002,
491210001,
501210001,
521210002,
551210001,
451220002,
471220001,
471220002,
491220001,
491220005,
511220005,
511220006,
551220007,
551220008,
481230003,
491230001,
501230001,
521230002,
551230005,
461240004,
471240001,
471240002,
471240003,
491240002,
501240016,
511240001,
511240002,
551240025,
471250007,
471250010,
481250001,
491250001,
501250001,
521250002,
541250002,
571250005,
461260003,
461260004,
461260005,
471260001,
471260004,
491260001,
511260001,
511260002,
481270006,
491270001,
491270009,
501270001,
521270002,
541270002,
561270002,
571270001,
581270001,
461280004,
511280001,
571280001,
471290001,
481290001,
481290004,
491290001,
491290010,
491290012,
491290013,
501290001,
501290017,
501290018,
501290025,
511290011,
511290012,
511290023,
521290001,
541290002,
551290010,
561290001,
571290002,
601290001,
601290003,
491300001,
491300002,
491300003,
501300002,
511300001,
531300001,
551300004,
561300030,
591300002,
491310001,
491310004,
501310001,
521310001,
521310033,
541310002,
561310002,
571310006,
581310001,
591310002,
501320006,
511320001,
521320006,
521320022,
531320003,
541320030,
571320004,
581320030,
491330001,
521330002,
531330016,
531330059,
531330065,
541330001,
561330002,
581330001,
591330003,
601330001,
611330005,
621330000,
511340002,
521340003,
531340005,
541340007,
551340003,
601340017,
611340000,
611340001,
521350010,
541350002,
551350010,
561350002,
581350004,
591350004,
601350001,
611350000,
611350003,
501360003,
531360006,
551360001,
561360005,
571360051,
611360000,
611360001,
631360001,
561370002,
581370002,
601370004,
501380003,
531380002,
551380003,
581380005,
591380005,
611380001,
571390021,
581390002,
601390002,
601390035,
611390001,
621390004,
631390003,
641390001,
511400003,
591400003,
591400015,
601400009,
611400008,
631400004,
601410002,
621410002,
631410001,
641410004,
651410001,
591420001,
591420024,
601420004,
611420012,
631420031,
641420019,
641420020,
651420003,
621430002,
621430043,
641430002,
651430001,
661430003,
551440004,
591440001,
651440004,
651440006,
651440007,
671440003,
641450002,
651450004,
661450002,
681450002,
571460001,
651460022,
651460026,
661460011,
691460001,
651470001,
661470002,
681470002,
691470001,
591480000,
591480001,
611480003,
651480001,
671480001,
671480012,
681480008,
651490001,
661490027,
671490001,
681490002,
631500001,
651500002,
671500001,
691500005,
581510001,
621510012,
631510002,
651510003,
671510001,
681510021,
691510001,
691510012,
701510001,
701510005,
701510010,
611520004,
611520014,
631520001,
631520016,
651520006,
671520001,
691520006,
691520018,
691520019,
701520006,
621530006,
641530003,
641530008,
651530003,
671530001,
691530001,
601540003,
611540001,
631540013,
651540001,
651540002,
711540015,
721540006,
641550006,
661550009,
671550002,
691550001,
711550001,
711550004,
611560002,
651560002,
651560004,
671560001,
671560012,
711560001,
721560004,
731560001,
641570012,
661570005,
711570001,
731570001,
731570004,
601580004,
611580002,
651580003,
651580023,
671580001,
671580009,
691580002,
691580003,
711580000,
731580001,
731580015,
621590006,
641590002,
661590009,
671590003,
731590001,
601600003,
671600001,
671600006,
691600002,
711600001,
611610005,
671610002,
681610014,
691610001,
711610004,
721610002,
671620003,
691620020,
711620008,
711620009,
751620001,
671630003,
751630001,
621640005,
671640003,
691640001,
771640001,
661650002,
751650001,
771650001,
641660009,
671660001,
691660006,
711660001,
711660002,
681670003,
711670001,
751670001,
671680001,
711680013,
771680001,
701690001,
711690001,
751690001,
771690001,
671700001,
711700008,
771700001,
711710001,
721710001,
771710001,
781710002,
691720011,
711720001,
711720005,
751720001,
771720001,
791720001,
771730000,
771730003,
791730001,
711740003,
771740001,
701750007,
711750053,
791750001,
701760005,
711760001,
731760012,
731760090,
791760001,
791760002,
691770000,
701770006,
711770029,
711770203,
721770048,
721770107,
791770002,
711780003,
721780005,
721780109,
731780000,
731780059,
731780094,
731780139,
711790006,
721790005,
721790046,
731790117,
741790002,
751790137,
791790007,
811790001,
721800006,
731800002,
721810025,
721810078,
761810001,
811810002,
721820009,
721820026,
731820001,
731820029,
741820062,
751820001,
761820029,
721830007,
731830032,
741830007,
751830004,
751830005,
751830070,
751830074,
761830002,
781830001,
811830002,
811830005,
821830001,
721840005,
751840005,
771840007,
781840034,
791840003,
741850006,
781850002,
791850001,
801850004,
811850003,
751860004,
771860001,
811860000,
811860005,
831860001,
791870002,
801870001,
811870002,
821870001,
831870002,
751880007,
811880001,
731890001,
751890033,
751890034,
761890001,
771890006,
771890085,
781890004,
781890005,
791890003,
791890006,
791890200,
801890002,
811890001,
821890001,
821890014,
831890002,
831890003,
731900002,
741900006,
751900003,
761900032,
771900002,
771900037,
791900014,
811900000,
811900001,
811900006,
831900000,
831900001,
761910001,
771910003,
771910071,
791910004,
801910035,
811910002,
821910002,
831910002,
831910005,
831910028,
751920002,
751920003,
761920047,
761920112,
771920003,
771920015,
791920004,
791920015,
811920002,
811920008,
821920011,
821920014,
821920017,
821920020,
821920021,
831920001,
841920006,
851920000,
851920001,
771930002,
781930005,
791930004,
801930003,
811930002,
821930001,
831930001,
841930001,
851930001,
851930002,
751940001,
751940002,
751940003,
771940007,
771940012,
791940003,
791940008,
811940001,
831940001,
831940002,
851940000,
851940001,
761950002,
761950004,
771950002,
781950007,
791950004,
791950055,
801950003,
811950002,
821950002,
831950001,
841950002,
851950001,
861950001,
751960001,
771960004,
791960003,
791960054,
811960006,
831960002,
831960003,
841960015,
761970001,
771970002,
781970009,
791970004,
801970004,
811970002,
821970002,
831970001,
841970002,
851970001,
861970001,
761980006,
761980010,
771980001,
791980051,
811980007,
811980012,
831980001,
831980003,
851980001,
871980001,
781990008,
791990006,
801990007,
811990003,
821990003,
831990001,
841990002,
861990001,
792000011,
812000010,
832000001,
832000003,
852000001,
852000003,
802010013,
812010003,
822010004,
832010001,
842010003,
862010001,
872010001,
882010000,
782020003,
822020014,
852020001,
852020002,
872020001,
822030006,
822030053,
832030006,
842030005,
862030001,
882030001,
812040029,
822040021,
832040008,
832040038,
852040001,
872040001,
872040002,
802050008,
822050009,
842050010,
842050017,
882050001,
812060045,
832060016,
872060001,
872060002,
892060001,
812070002,
822070003,
832070036,
842070014,
862070007,
882070001,
802080004,
832080018,
802100002,
802100005,
832100002,
822110014,
832110021,
842110015,
852110076,
872110013,
872110019,
832120005,
832120012,
842120030,
852120004,
882130005,
822140004,
852140006,
862140004,
862140005,
872140001,
902140004,
832150009,
862150013,
902150003,
822160004,
872160001,
922160001,
832170005,
892170010,
902170001,
912170001,
872180002,
922180001,
892220001,
912340002,
922350001,
922350174,
942350010,
932360001,
952360001,
942370003,
922380119,
932380128,
942380047,
942380052,
952380001,
942390106,
942390111,
952390011,
932400001,
942400102,
952400057,
962400002,
962400003,
942410106,
942410107,
952410075,
962410007,
932420007,
942420044,
942420045,
952420002,
952420141,
962420004,
962420005,
972420002,
972420003,
942440032,
952440001,
952440112,
952440113,
962440009,
962440013,
962440014,
972440004,
982440002,
942450024,
952450021,
962450061,
972450003,
1012450001,
952460001,
952460008,
972460000,
982460002,
992460000,
1012460000,
1012460001,
1002470001,
1012470001,
972480001,
1002480006,
992500001,
1002500001,
1002500002,
1022500001,
1022510002,
1002530008,
1022530003,
1022530030,
1022530031,
1022530032,
1032530000,
1032530001,
992540002,
1012540000,
1012540001,
1022540011,
1032550001,
1032550027,
992560001,
1002560026,
1042560007,
1042560009,
1042560012,
1042570002,
1052570002,
1012580001,
1052580001,
1042610001,
1072620001,
1062630003,
1062650001,
1082650001,
1082670002,
1102700001,
1102710001,
1082770001,
};
int map_metastable [TOTAL_STATE_MAPS] = {1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
2,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
2,
3,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
3,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
3,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
2,
1,
2,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
3,
1,
1,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
3,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
2,
3,
4,
1,
2,
3,
4,
1,
2,
3,
1,
1,
1,
1,
1,
1,
2,
1,
2,
3,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
3,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
2,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
2,
3,
1,
2,
1,
2,
1,
1,
1,
2,
3,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
2,
1,
2,
1,
1,
1,
2,
1,
1,
1,
2,
1,
1,
1,
2,
1,
2,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
2,
1,
1,
1,
2,
1,
2,
3,
4,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
3,
1,
1,
1,
2,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
3,
4,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
2,
1,
2,
1,
2,
3,
1,
1,
1,
2,
1,
2,
1,
1,
1,
1,
1,
2,
1,
1,
2,
3,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
2,
3,
1,
2,
1,
2,
1,
2,
1,
2,
1,
2,
1,
2,
3,
4,
5,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
2,
3,
1,
2,
1,
2,
1,
1,
2,
1,
2,
1,
2,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
2,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
2,
1,
1,
1,
2,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
2,
1,
2,
2,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
2,
1,
1,
1,
1,
1,
2,
1,
2,
1,
1,
2,
1,
2,
1,
2,
1,
2,
1,
2,
1,
1,
2,
3,
1,
2,
3,
1,
1,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
2,
1,
1,
1,
1,
1,
1,
2,
1,
1,
1,
1,
2,
3,
4,
1,
1,
1,
1,
2,
1,
1,
2,
1,
1,
1,
2,
3,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
};
}
}
#endif
//
// end of /home/mouginot/work/app/pyne/src/state_map.cpp
//


//
// start of /home/mouginot/work/app/pyne/src/nucname.cpp
//
// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// name is for letters  as well (U-235).
// MCNP is for numerals without the meta-stable flag (92235), as used in MCNP.

#ifndef PYNE_IS_AMALGAMATED
#include "nucname.h"
#include "state_map.cpp"
#endif


/*** Constructs the LL to zz Dictionary ***/
pyne::nucname::name_zz_t pyne::nucname::get_name_zz() {
  pyne::nucname::name_zz_t lzd;

  lzd["Be"] = 04;
  lzd["Ba"] = 56;
  lzd["Bh"] = 107;
  lzd["Bi"] = 83;
  lzd["Bk"] = 97;
  lzd["Br"] = 35;
  lzd["Ru"] = 44;
  lzd["Re"] = 75;
  lzd["Rf"] = 104;
  lzd["Rg"] = 111;
  lzd["Ra"] = 88;
  lzd["Rb"] = 37;
  lzd["Rn"] = 86;
  lzd["Rh"] = 45;
  lzd["Tm"] = 69;
  lzd["H"] = 01;
  lzd["P"] = 15;
  lzd["Ge"] = 32;
  lzd["Gd"] = 64;
  lzd["Ga"] = 31;
  lzd["Os"] = 76;
  lzd["Hs"] = 108;
  lzd["Zn"] = 30;
  lzd["Ho"] = 67;
  lzd["Hf"] = 72;
  lzd["Hg"] = 80;
  lzd["He"] = 02;
  lzd["Pr"] = 59;
  lzd["Pt"] = 78;
  lzd["Pu"] = 94;
  lzd["Pb"] = 82;
  lzd["Pa"] = 91;
  lzd["Pd"] = 46;
  lzd["Po"] = 84;
  lzd["Pm"] = 61;
  lzd["C"] = 6;
  lzd["K"] = 19;
  lzd["O"] = 8;
  lzd["S"] = 16;
  lzd["W"] = 74;
  lzd["Eu"] = 63;
  lzd["Es"] = 99;
  lzd["Er"] = 68;
  lzd["Md"] = 101;
  lzd["Mg"] = 12;
  lzd["Mo"] = 42;
  lzd["Mn"] = 25;
  lzd["Mt"] = 109;
  lzd["U"] = 92;
  lzd["Fr"] = 87;
  lzd["Fe"] = 26;
  lzd["Fm"] = 100;
  lzd["Ni"] = 28;
  lzd["No"] = 102;
  lzd["Na"] = 11;
  lzd["Nb"] = 41;
  lzd["Nd"] = 60;
  lzd["Ne"] = 10;
  lzd["Zr"] = 40;
  lzd["Np"] = 93;
  lzd["B"] = 05;
  lzd["Co"] = 27;
  lzd["Cm"] = 96;
  lzd["F"] = 9;
  lzd["Ca"] = 20;
  lzd["Cf"] = 98;
  lzd["Ce"] = 58;
  lzd["Cd"] = 48;
  lzd["V"] = 23;
  lzd["Cs"] = 55;
  lzd["Cr"] = 24;
  lzd["Cu"] = 29;
  lzd["Sr"] = 38;
  lzd["Kr"] = 36;
  lzd["Si"] = 14;
  lzd["Sn"] = 50;
  lzd["Sm"] = 62;
  lzd["Sc"] = 21;
  lzd["Sb"] = 51;
  lzd["Sg"] = 106;
  lzd["Se"] = 34;
  lzd["Yb"] = 70;
  lzd["Db"] = 105;
  lzd["Dy"] = 66;
  lzd["Ds"] = 110;
  lzd["La"] = 57;
  lzd["Cl"] = 17;
  lzd["Li"] = 03;
  lzd["Tl"] = 81;
  lzd["Lu"] = 71;
  lzd["Lr"] = 103;
  lzd["Th"] = 90;
  lzd["Ti"] = 22;
  lzd["Te"] = 52;
  lzd["Tb"] = 65;
  lzd["Tc"] = 43;
  lzd["Ta"] = 73;
  lzd["Ac"] = 89;
  lzd["Ag"] = 47;
  lzd["I"] = 53;
  lzd["Ir"] = 77;
  lzd["Am"] = 95;
  lzd["Al"] = 13;
  lzd["As"] = 33;
  lzd["Ar"] = 18;
  lzd["Au"] = 79;
  lzd["At"] = 85;
  lzd["In"] = 49;
  lzd["Y"] = 39;
  lzd["N"] = 07;
  lzd["Xe"] = 54;
  lzd["Cn"] = 112;
  lzd["Fl"] = 114;
  lzd["Lv"] = 116;

  return lzd;
}
pyne::nucname::name_zz_t pyne::nucname::name_zz = pyne::nucname::get_name_zz();


/*** Constructs zz to LL dictionary **/
pyne::nucname::zzname_t pyne::nucname::get_zz_name()
{
  zzname_t zld;
  for (name_zz_iter i = name_zz.begin(); i != name_zz.end(); i++)
  {
    zld[i->second] = i->first;
  }
  return zld;
}
pyne::nucname::zzname_t pyne::nucname::zz_name = pyne::nucname::get_zz_name();



/*** Constructs the fluka to zz Dictionary ***/
pyne::nucname::name_zz_t pyne::nucname::get_fluka_zz() {
  pyne::nucname::name_zz_t fzd;

  fzd["BERYLLIU"] = 40000000;
  fzd["BARIUM"]   = 560000000;
  fzd["BOHRIUM"]  = 1070000000;   // No fluka
  fzd["BISMUTH"]  = 830000000;
  fzd["BERKELIU"] = 970000000;    // No fluka
  fzd["BROMINE"]  = 350000000;
  fzd["RUTHENIU"] = 440000000;    // No fluka
  fzd["RHENIUM"]  = 750000000;
  fzd["RUTHERFO"] = 1040000000;
  fzd["ROENTGEN"] = 1110000000;
  fzd["RADIUM"]   = 880000000;    // No fluka
  fzd["RUBIDIUM"] = 370000000;    // No fluka
  fzd["RADON"]    = 860000000;    // no fluka
  fzd["RHODIUM"]  = 450000000;    // no fluka
  fzd["THULIUM"]  = 690000000;    // no fluka
  fzd["HYDROGEN"] = 10000000;
  fzd["PHOSPHO"]  = 150000000;
  fzd["GERMANIU"] = 320000000;
  fzd["GADOLINI"] = 640000000;
  fzd["GALLIUM"]  = 310000000;
  fzd["OSMIUM"]   = 760000000;    // no fluka
  fzd["HASSIUM"]  = 1080000000;
  fzd["ZINC"]     = 300000000;
  fzd["HOLMIUM"]  = 670000000;    // no fluka
  fzd["HAFNIUM"]  = 720000000;
  fzd["MERCURY"]  = 800000000;
  fzd["HELIUM"]   = 20000000;
  fzd["PRASEODY"] = 590000000;   // no fluka
  fzd["PLATINUM"] = 780000000;
  fzd["239-PU"]   = 940000000;   // "239-PU"
  fzd["LEAD"]     = 820000000;
  fzd["PROTACTI"] = 910000000;   // no fluka
  fzd["PALLADIU"] = 460000000;   // no fluka
  fzd["POLONIUM"] = 840000000;   // no fluka
  fzd["PROMETHI"] = 610000000;   // no fluka
  fzd["CARBON"]   = 60000000;
  fzd["POTASSIU"] = 190000000;
  fzd["OXYGEN"]   = 80000000;
  fzd["SULFUR"]   = 160000000;
  fzd["TUNGSTEN"] = 740000000;
  fzd["EUROPIUM"] = 630000000;
  fzd["EINSTEIN"] = 990000000;   // no fluka
  fzd["ERBIUM"]   = 680000000;   // no fluka
  fzd["MENDELEV"] = 1010000000;  // no fluka
  fzd["MAGNESIU"] = 120000000;
  fzd["MOLYBDEN"] = 420000000;
  fzd["MANGANES"] = 250000000;
  fzd["MEITNERI"] = 1090000000;  // no fluka
  fzd["URANIUM"]  = 920000000;
  fzd["FRANCIUM"] = 870000000;   // no fluka
  fzd["IRON"]     = 260000000;
  fzd["FERMIUM"]  = 1000000000;  // no fluka
  fzd["NICKEL"]   = 280000000;
  fzd["NITROGEN"] = 70000000;
  fzd["NOBELIUM"] = 1020000000;  // no fluka
  fzd["SODIUM"]   = 110000000;
  fzd["NIOBIUM"]  = 410000000;
  fzd["NEODYMIU"] = 600000000;
  fzd["NEON"]     = 100000000;
  fzd["ZIRCONIU"] = 400000000;
  fzd["NEPTUNIU"] = 930000000;   // no fluka
  fzd["BORON"]    = 50000000;
  fzd["COBALT"]   = 270000000;
  fzd["CURIUM"]   = 960000000;   // no fluka
  fzd["FLUORINE"] = 90000000;
  fzd["CALCIUM"]  = 200000000;
  fzd["CALIFORN"] = 980000000;   // no fluka
  fzd["CERIUM"]   = 580000000;
  fzd["CADMIUM"]  = 480000000;
  fzd["VANADIUM"] = 230000000;
  fzd["CESIUM"]   = 550000000;
  fzd["CHROMIUM"] = 240000000;
  fzd["COPPER"]   = 290000000;
  fzd["STRONTIU"] = 380000000;
  fzd["KRYPTON"]  = 360000000;
  fzd["SILICON"]  = 140000000;
  fzd["TIN"]      = 500000000;
  fzd["SAMARIUM"] = 620000000;
  fzd["SCANDIUM"] = 210000000;
  fzd["ANTIMONY"] = 510000000;
  fzd["SEABORGI"] = 1060000000;  // no fluka
  fzd["SELENIUM"] = 340000000;   // no fluka
  fzd["YTTERBIU"] = 700000000;   // no fluka
  fzd["DUBNIUM"]  = 1050000000;  // no fluka
  fzd["DYSPROSI"] = 660000000;   // no fluka
  fzd["DARMSTAD"] = 1100000000;  // no fluka
  fzd["LANTHANU"] = 570000000;
  fzd["CHLORINE"] = 170000000;
  fzd["LITHIUM"]  = 030000000;
  fzd["THALLIUM"] = 810000000;   // no fluka
  fzd["LUTETIUM"] = 710000000;   // no fluka
  fzd["LAWRENCI"] = 1030000000;  // no fluka
  fzd["THORIUM"]  = 900000000;   // no fluka
  fzd["TITANIUM"] = 220000000;
  fzd["TELLURIU"] = 520000000;   // no fluka
  fzd["TERBIUM"]  = 650000000;
  fzd["99-TC"]    = 430000000;   // "99-TC"
  fzd["TANTALUM"] = 730000000;
  fzd["ACTINIUM"] = 890000000;   // no fluka
  fzd["SILVER"]   = 470000000;
  fzd["IODINE"]   = 530000000;
  fzd["IRIDIUM"]  = 770000000;
  fzd["241-AM"]   = 950000000;   // "241-AM"
  fzd["ALUMINUM"] = 130000000;
  fzd["ARSENIC"]  = 330000000;
  fzd["ARGON"]    = 180000000;
  fzd["GOLD"]     = 790000000;
  fzd["ASTATINE"] = 850000000;   // no fluka
  fzd["INDIUM"]   = 490000000;
  fzd["YTTRIUM"]  = 390000000;
  fzd["XENON"]    = 540000000;
  fzd["COPERNIC"] = 1120000000;  // no fluka
  fzd["UNUNQUAD"] = 1140000000;  // no fluka:  UNUNQUADIUM,  "Flerovium"
  fzd["UNUNHEXI"] = 1160000000;  // no fluka:  UNUNHEXIUM , "Livermorium"
  fzd["HYDROG-1"] = 10010000;
  fzd["DEUTERIU"] = 10020000;
  fzd["TRITIUM"]  = 10040000;
  fzd["HELIUM-3"] = 20030000;
  fzd["HELIUM-4"] = 20040000;
  fzd["LITHIU-6"] = 30060000;
  fzd["LITHIU-7"] = 30070000;
  fzd["BORON-10"] = 50100000;
  fzd["BORON-11"] = 50110000;
  fzd["90-SR"]    = 380900000;   // fluka "90-SR"
  fzd["129-I"]    = 531290000;   // fluka "129-I"
  fzd["124-XE"]   = 541240000;   // fluka "124-XE"
  fzd["126-XE"]   = 541260000;   // fluka "126-XE"
  fzd["128-XE"]   = 541280000;   // fluka "128-XE"
  fzd["130-XE"]   = 541300000;   // fluka "130-XE"
  fzd["131-XE"]   = 541310000;   // fluka "131-XE"
  fzd["132-XE"]   = 541320000;   // fluka "132-XE"
  fzd["134-XE"]   = 541340000;   // fluka "134-XE"
  fzd["135-XE"]   = 541350000;   // fluka "135-XE"
  fzd["136-XE"]   = 541360000;   // fluka "136-XE"
  fzd["135-CS"]   = 551350000;   // fluka "135-CS"
  fzd["137-CS"]   = 551370000;   // fluka "137-CS"
  fzd["230-TH"]   = 902300000;   // fluka "230-TH"
  fzd["232-TH"]   = 902320000;   // fluka "232-TH"
  fzd["233-U"]    = 922330000;   // fluka "233-U"
  fzd["234-U"]    = 922340000;   // fluka "234-U"
  fzd["235-U"]    = 922350000;   // fluka "235-U"
  fzd["238-U"]    = 922380000;   // fluka "238-U"

  return fzd;
}
pyne::nucname::name_zz_t pyne::nucname::fluka_zz = pyne::nucname::get_fluka_zz();


/*** Constructs zz to fluka dictionary **/
pyne::nucname::zzname_t pyne::nucname::get_zz_fluka()
{
  zzname_t zfd;
  for (name_zz_iter i = fluka_zz.begin(); i != fluka_zz.end(); i++)
  {
    zfd[i->second] = i->first;
  }
  return zfd;
}
pyne::nucname::zzname_t pyne::nucname::zz_fluka = pyne::nucname::get_zz_fluka();



/******************************************/
/*** Define useful elemental group sets ***/
/******************************************/

pyne::nucname::zz_group pyne::nucname::name_to_zz_group(pyne::nucname::name_group eg)
{
  zz_group zg;
  for (name_group_iter i = eg.begin(); i != eg.end(); i++)
    zg.insert(name_zz[*i]);
  return zg;
}

// Lanthanides
pyne::nucname::name_t pyne::nucname::LAN_array[15] = {"La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"};
pyne::nucname::name_group pyne::nucname::LAN (pyne::nucname::LAN_array,
                                              pyne::nucname::LAN_array+15);
pyne::nucname::zz_group pyne::nucname::lan = \
  pyne::nucname::name_to_zz_group(pyne::nucname::LAN);

// Actinides
pyne::nucname::name_t pyne::nucname::ACT_array[15] = {"Ac", "Th", "Pa", "U",
  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
pyne::nucname::name_group pyne::nucname::ACT (pyne::nucname::ACT_array, pyne::nucname::ACT_array+15);
pyne::nucname::zz_group pyne::nucname::act = pyne::nucname::name_to_zz_group(pyne::nucname::ACT);

// Transuarnics
pyne::nucname::name_t pyne::nucname::TRU_array[22] = {"Np", "Pu", "Am", "Cm",
  "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
  "Ds", "Rg", "Cn", "Fl", "Lv"};
pyne::nucname::name_group pyne::nucname::TRU (pyne::nucname::TRU_array,
                                              pyne::nucname::TRU_array+22);
pyne::nucname::zz_group pyne::nucname::tru = \
  pyne::nucname::name_to_zz_group(pyne::nucname::TRU);

//Minor Actinides
pyne::nucname::name_t pyne::nucname::MA_array[10] = {"Np", "Am", "Cm", "Bk",
  "Cf", "Es", "Fm", "Md", "No", "Lr"};
pyne::nucname::name_group pyne::nucname::MA (pyne::nucname::MA_array,
                                             pyne::nucname::MA_array+10);
pyne::nucname::zz_group pyne::nucname::ma = \
  pyne::nucname::name_to_zz_group(pyne::nucname::MA);

//Fission Products
pyne::nucname::name_t pyne::nucname::FP_array[88] = {"Ag", "Al", "Ar", "As",
  "At", "Au", "B",  "Ba", "Be", "Bi", "Br", "C",  "Ca", "Cd", "Ce", "Cl", "Co",
  "Cr", "Cs", "Cu", "Dy", "Er", "Eu", "F",  "Fe", "Fr", "Ga", "Gd", "Ge", "H",
  "He", "Hf", "Hg", "Ho", "I",  "In", "Ir", "K",  "Kr", "La", "Li", "Lu", "Mg",
  "Mn", "Mo", "N",  "Na", "Nb", "Nd", "Ne", "Ni", "O",  "Os", "P",  "Pb", "Pd",
  "Pm", "Po", "Pr", "Pt", "Ra", "Rb", "Re", "Rh", "Rn", "Ru", "S",  "Sb", "Sc",
  "Se", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Ti", "Tl", "Tm", "V",
  "W",  "Xe", "Y",  "Yb", "Zn", "Zr"};
pyne::nucname::name_group pyne::nucname::FP (pyne::nucname::FP_array,
                                             pyne::nucname::FP_array+88);
pyne::nucname::zz_group pyne::nucname::fp = \
  pyne::nucname::name_to_zz_group(pyne::nucname::FP);


/***************************/
/*** isnuclide functions ***/
/***************************/

bool pyne::nucname::isnuclide(std::string nuc) {
  int n;
  try {
    n = id(nuc);
  }
  catch(NotANuclide) {
    return false;
  }
  catch(IndeterminateNuclideForm) {
    return false;
  }
  return isnuclide(n);
}

bool pyne::nucname::isnuclide(const char * nuc) {
  return isnuclide(std::string(nuc));
}

bool pyne::nucname::isnuclide(int nuc) {
  int n;
  try {
    n = id(nuc);
  }
  catch(NotANuclide) {
    return false;
  }
  catch(IndeterminateNuclideForm) {
    return false;
  }
  if (n <= 10000000)
    return false;
  int zzz = n / 10000000;
  int aaa = (n % 10000000) / 10000;
  if (aaa == 0)
    return false;  // is element
  else if (aaa < zzz)
    return false;
  return true;
}



/********************/
/*** id functions ***/
/********************/
int pyne::nucname::id(int nuc) {
  if (nuc < 0)
    throw NotANuclide(nuc, "");

  int newnuc;
  int zzz = nuc / 10000000;     // ZZZ ?
  int aaassss = nuc % 10000000; // AAA-SSSS ?
  int aaa = aaassss / 10000;    // AAA ?
  int ssss = aaassss % 10000;   // SSSS ?
  // Nuclide must already be in id form
  if (0 < zzz && zzz <= aaa && aaa <= zzz * 7) {
    // Normal nuclide
    if (5 < ssss){
    // Unphysical metastable state warning
     warning("You have indicated a metastable state of " + pyne::to_str(ssss) + ". Metastable state above 5, possibly unphysical. ");
    }
    return nuc;
  } else if (aaassss == 0 && 0 < zz_name.count(zzz)) {
    // Natural elemental nuclide:  ie for Uranium = 920000000
    return nuc;
  } else if (nuc < 1000 && 0 < zz_name.count(nuc))
    //  Gave Z-number
    return nuc * 10000000;

  // Not in id form, try  ZZZAAAM form.
  zzz = nuc / 10000;     // ZZZ ?
  aaassss = nuc % 10000; // AAA-SSSS ?
  aaa = aaassss / 10;    // AAA ?
  ssss = nuc % 10;       // SSSS ?
  if (zzz <= aaa && aaa <= zzz * 7) {
    // ZZZAAAM nuclide
    if (5 < ssss){
    // Unphysical metastable state warning
      warning("You have indicated a metastable state of " + pyne::to_str(ssss) + ". Metastable state above 5, possibly unphysical. ");
    }
    return (zzz*10000000) + (aaa*10000) + (nuc%10);
  } else if (aaa <= zzz && zzz <= aaa * 7 && 0 < zz_name.count(aaa)) {
    // Cinder-form (aaazzzm), ie 2350920
    if (5 < ssss){
    // Unphysical metastable state warning
      warning("You have indicated a metastable state of " + pyne::to_str(ssss) + ". Metastable state above 5, possibly unphysical. ");
    }
    return (aaa*10000000) + (zzz*10000) + (nuc%10);
  }
  //else if (aaassss == 0 && 0 == zz_name.count(nuc/1000) && 0 < zz_name.count(zzz))
  else if (aaassss == 0 && 0 < zz_name.count(zzz)) {
    // zzaaam form natural nuclide
    return zzz * 10000000;
  }

  if (nuc >= 1000000){
    // From now we assume no metastable info has been given.
    throw IndeterminateNuclideForm(nuc, "");
  }

  // Nuclide is not in zzaaam form,
  // Try MCNP form, ie zzaaa
  // This is the same form as SZA for the 0th state.
  zzz = nuc / 1000;
  aaa = nuc % 1000;
  if (zzz <= aaa) {
    if (aaa - 400 < 0) {
      if (nuc == 95242)
        return nuc * 10000 + 1;  // special case MCNP Am-242m
      else
        return nuc * 10000;  // Nuclide in normal MCNP form
    } else {
      // Nuclide in MCNP metastable form
      if (nuc == 95642)
        return (95642 - 400)*10000;  // special case MCNP Am-242
      nuc = ((nuc - 400) * 10000) + 1;
      while (3.0 < (float ((nuc/10000)%1000) / float (nuc/10000000)))
        nuc -= 999999;
      return nuc;
    }
  } else if (aaa == 0 && 0 < zz_name.count(zzz)) {
    // MCNP form natural nuclide
    return zzz * 10000000;
  }

  // Not a normal nuclide, might be a
  // Natural elemental nuclide.
  // ie 92 for Uranium = 920000
  if (0 < zz_name.count(nuc))
    return nuc * 10000000;
  throw IndeterminateNuclideForm(nuc, "");
}

int pyne::nucname::id(const char * nuc) {
  std::string newnuc (nuc);
  return id(newnuc);
}

int pyne::nucname::id(std::string nuc) {
  size_t npos = std::string::npos;
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int newnuc;
  std::string elem_name;
  int dash1 = nuc.find("-");
  int dash2;
  if (dash1 == npos)
    dash2 = npos;
  else
    dash2 = nuc.find("-", dash1+1);

  // nuc must be at least 4 characters or greater if it is in ZZLLAAAM form.
  if (nuc.length() >= 5 && dash1 != npos && dash2 != npos) {
    // Nuclide most likely in ZZLLAAAM Form, only form that contains two "-"'s.
    std::string zz = nuc.substr(0, dash1);
    std::string ll = nuc.substr(dash1+1, dash2);
    int zz_int = to_int(zz);
    // Verifying that the LL and ZZ point to the same element as secondary
    if(znum(ll) != zz_int)
      throw NotANuclide(nuc, "mismatched znum and chemical symbol");
    return zzllaaam_to_id(nuc);
  }

  // Get the string into a regular form
  std::string nucstr = pyne::to_upper(nuc);
  nucstr = pyne::remove_substring(nucstr, "-");
  int nuclen = nucstr.length();

  if (pyne::contains_substring(pyne::digits, nucstr.substr(0, 1))) {
    if (pyne::contains_substring(pyne::digits, nucstr.substr(nuclen-1, nuclen))) {
      // Nuclide must actually be an integer that
      // just happens to be living in string form.
      newnuc = pyne::to_int(nucstr);
      newnuc = id(newnuc);
    } else {
      // probably in NIST-like form (242Am)
      // Here we know we have both digits and letters
      std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);
      newnuc = pyne::to_int(anum_str) * 10000;

      // Add the Z-number
      elem_name = pyne::remove_characters(nucstr, pyne::digits);
      elem_name = pyne::capitalize(elem_name);
      if (0 < name_zz.count(elem_name))
        newnuc = (10000000 * name_zz[elem_name]) + newnuc;
      else
        throw NotANuclide(nucstr, newnuc);
    }
  } else if (pyne::contains_substring(pyne::alphabet, nucstr.substr(0, 1))) {
    // Nuclide is probably in name form, or some variation therein
    std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);

    // natural element form, a la 'U' -> 920000000
    if (anum_str.empty()) {
      elem_name = pyne::capitalize(nucstr);
      if (0 < name_zz.count(elem_name))
        return 10000000 * name_zz[elem_name];
    }

    int anum = pyne::to_int(anum_str);

    // bad form
    if (anum < 0)
      throw NotANuclide(nucstr, anum);

    // Figure out if we are meta-stable or not
    std::string end_char = pyne::last_char(nucstr);
    if (end_char == "M")
      newnuc = (10000 * anum) + 1;
    else if (pyne::contains_substring(pyne::digits, end_char))
      newnuc = (10000 * anum);
    else
      throw NotANuclide(nucstr, newnuc);

    // Add the Z-number
    elem_name = pyne::remove_characters(nucstr.substr(0, nuclen-1), pyne::digits);
    elem_name = pyne::capitalize(elem_name);
    if (0 < name_zz.count(elem_name))
      newnuc = (10000000 * name_zz[elem_name]) + newnuc;
    else
      throw NotANuclide(nucstr, newnuc);
  } else {
    // Clearly not a nuclide
    throw NotANuclide(nuc, nucstr);
  }
  return newnuc;
}


/***************************/
/*** iselement functions ***/
/***************************/

bool pyne::nucname::iselement(std::string nuc) {
  int n;
  try {
    n = id(nuc);
  }
  catch(NotANuclide) {
    return false;
  }
  return iselement(n);
}

bool pyne::nucname::iselement(const char * nuc) {
  return iselement(std::string(nuc));
}

bool pyne::nucname::iselement(int nuc) {
  int n;
  try {
    n = id(nuc);
  }
  catch(NotANuclide) {
    return false;
  }

  if (n < 10000000)
    return false;
  int zzz = znum(n);
  int aaa = anum(n);
  if (zzz > 0 && aaa == 0)
    return true;  // is element
  return false;
}

/**********************/
/*** name functions ***/
/**********************/
std::string pyne::nucname::name(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add LL
  newnuc += zz_name[zzz];

  // Add A-number
  if (0 < aaa)
    newnuc += pyne::to_str(aaa);

  // Add meta-stable flag
  if (0 < ssss)
    newnuc += "M";

  return newnuc;
}



std::string pyne::nucname::name(const char * nuc) {
  std::string newnuc (nuc);
  return name(newnuc);
}


std::string pyne::nucname::name(std::string nuc) {
  return name(id(nuc));
}


/**********************/
/*** znum functions ***/
/**********************/
int pyne::nucname::znum(int nuc) {
  return id(nuc) / 10000000;
}

int pyne::nucname::znum(const char * nuc) {
  return id(nuc) / 10000000;
}

int pyne::nucname::znum(std::string nuc) {
  return id(nuc) / 10000000;
}

/**********************/
/*** anum functions ***/
/**********************/
int pyne::nucname::anum(int nuc) {
  return (id(nuc) / 10000) % 1000;
}

int pyne::nucname::anum(const char * nuc) {
  return (id(nuc) / 10000) % 1000;
}

int pyne::nucname::anum(std::string nuc) {
  return (id(nuc) / 10000) % 1000;
}

/**********************/
/*** snum functions ***/
/**********************/
int pyne::nucname::snum(int nuc) {
  return id(nuc) % 10000;
}

int pyne::nucname::snum(const char * nuc) {
  return id(nuc) % 10000;
}

int pyne::nucname::snum(std::string nuc) {
  return id(nuc) % 10000;
}

/************************/
/*** zzaaam functions ***/
/************************/
int pyne::nucname::zzaaam(int nuc) {
  int nucid = id(nuc);
  int zzzaaa = nucid / 10000;
  int ssss = nucid % 10000;
  if (10 <= ssss)
    ssss = 9;
  return zzzaaa*10 + ssss;
}


int pyne::nucname::zzaaam(const char * nuc) {
  std::string newnuc (nuc);
  return zzaaam(newnuc);
}


int pyne::nucname::zzaaam(std::string nuc) {
  return zzaaam(id(nuc));
}


int pyne::nucname::zzaaam_to_id(int nuc) {
  return (nuc/10)*10000 + (nuc%10);
}


int pyne::nucname::zzaaam_to_id(const char * nuc) {
  return zzaaam_to_id(std::string(nuc));
}


int pyne::nucname::zzaaam_to_id(std::string nuc) {
  return zzaaam_to_id(pyne::to_int(nuc));
}

/************************/
/*** zzzaaa functions ***/
/************************/
int pyne::nucname::zzzaaa(int nuc) {
  int nucid = id(nuc);
  int zzzaaa = nucid/10000;

  return zzzaaa;
}


int pyne::nucname::zzzaaa(const char * nuc) {
  std::string newnuc (nuc);
  return zzzaaa(newnuc);
}


int pyne::nucname::zzzaaa(std::string nuc) {
  return zzzaaa(id(nuc));
}


int pyne::nucname::zzzaaa_to_id(int nuc) {
  return (nuc)*10000;
}


int pyne::nucname::zzzaaa_to_id(const char * nuc) {
  return zzzaaa_to_id(std::string(nuc));
}


int pyne::nucname::zzzaaa_to_id(std::string nuc) {
  return zzzaaa_to_id(pyne::to_int(nuc));
}

/*************************/
/*** zzllaaam functions ***/
/*************************/
std::string pyne::nucname::zzllaaam(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";

  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int zzz = nucid / 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);
  //Adding ZZ
  newnuc += pyne::to_str(zzz);
  newnuc += "-";
  // Add LL
  newnuc += zz_name[zzz];
  // Add required dash
  newnuc += "-";
  // Add AAA
  if (0 < aaassss)
    newnuc += pyne::to_str(aaa);
  // Add meta-stable flag
  if (0 < ssss)
    newnuc += "m";
  return newnuc;
}


std::string pyne::nucname::zzllaaam(const char * nuc) {
  std::string newnuc (nuc);
  return zzllaaam(newnuc);
}


std::string pyne::nucname::zzllaaam(std::string nuc) {
  return zzllaaam(id(nuc));
}


int pyne::nucname::zzllaaam_to_id(const char * nuc) {
  return zzllaaam_to_id(std::string(nuc));
}


int pyne::nucname::zzllaaam_to_id(std::string nuc) {
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int nucid;
  std::string elem_name;

  // Get the string into a regular form
  std::string nucstr = pyne::to_upper(nuc);
  // Removing first two characters (redundant), for 1 digit nuclides, such
  // as 2-He-4, the first slash will be removed, and the second attempt to
  // remove the second slash will do nothing.
  nucstr.erase(0,2);
  nucstr = pyne::remove_substring(nucstr, "-");
  // Does nothing if nuclide is short, otherwise removes the second "-" instance
  nucstr = pyne::remove_substring(nucstr, "-");
  int nuclen = nucstr.length();

  // Nuclide is probably in name form, or some variation therein
  std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);

  // natural element form, a la 'U' -> 920000000
  if (anum_str.empty() || pyne::contains_substring(nucstr, "NAT")) {
    elem_name = pyne::capitalize(pyne::remove_substring(nucstr, "NAT"));
    if (0 < name_zz.count(elem_name))
      return 10000000 * name_zz[elem_name];
  }
  int anum = pyne::to_int(anum_str);

  // Figure out if we are meta-stable or not
  std::string end_char = pyne::last_char(nucstr);
  if (end_char == "M")
    nucid = (10000 * anum) + 1;
  else if (pyne::contains_substring(pyne::digits, end_char))
    nucid = (10000 * anum);
  else
    throw NotANuclide(nucstr, nucid);

  // Add the Z-number
  elem_name = pyne::remove_characters(nucstr.substr(0, nuclen-1), pyne::digits);
  elem_name = pyne::capitalize(elem_name);
  if (0 < name_zz.count(elem_name))
    nucid = (10000000 * name_zz[elem_name]) + nucid;
  else
    throw NotANuclide(nucstr, nucid);
  return nucid;
}

/**********************/
/*** mcnp functions ***/
/**********************/
int pyne::nucname::mcnp(int nuc) {
  nuc = id(nuc);
  int ssss = nuc % 10000;
  int newnuc = nuc / 10000;

  // special case Am242(m)
  if (newnuc == 95242 && ssss < 2)
    ssss = (ssss + 1) % 2;

  // Handle the crazy MCNP meta-stable format
  if (0 != ssss && ssss < 10)
    newnuc += 300 + (ssss * 100);

  return newnuc;
}



int pyne::nucname::mcnp(const char * nuc) {
  std::string newnuc (nuc);
  return mcnp(newnuc);
}



int pyne::nucname::mcnp(std::string nuc) {
  return mcnp(id(nuc));
}

//
// MCNP -> id
//
int pyne::nucname::mcnp_to_id(int nuc) {
  int zzz = nuc / 1000;
  int aaa = nuc % 1000;
  if (zzz == 0)
    throw NotANuclide(nuc, "not in the MCNP format");
  else if (zzz <= aaa) {
    if (aaa - 400 < 0) {
      if (nuc == 95242)
        return nuc * 10000 + 1;  // special case MCNP Am-242m
      else
        return nuc * 10000;  // Nuclide in normal MCNP form
    } else {
      // Nuclide in MCNP metastable form
      if (nuc == 95642)
        return (95642 - 400)*10000;  // special case MCNP Am-242
      nuc = ((nuc - 400) * 10000) + 1;
      while (3.0 < (float ((nuc/10000)%1000) / float (nuc/10000000)))
        nuc -= 999999;
      return nuc;
    }
  } else if (aaa == 0)
    // MCNP form natural nuclide
    return zzz * 10000000;
  throw IndeterminateNuclideForm(nuc, "");
}


int pyne::nucname::mcnp_to_id(const char * nuc) {
  return mcnp_to_id(std::string(nuc));
}


int pyne::nucname::mcnp_to_id(std::string nuc) {
  return mcnp_to_id(pyne::to_int(nuc));
}

/************************/
/*** openmc functions ***/
/************************/
std::string pyne::nucname::openmc(int nuc) {
  std::string nucname = name(nuc);

  // check aaa value
  if (iselement(nuc)) {
    nucname.append("0");
  }

  // format metadata
  if ('M' == nucname.back()) {
    nucname.back() = '_';
    nucname.append("m");
    int meta_id = snum(nuc);
    std::string meta_str = std::to_string(meta_id);
    nucname.append(meta_str);
  }
  return nucname;
}

std::string pyne::nucname::openmc(const char * nuc) {
  std::string newnuc (nuc);
  return openmc(newnuc);
}

std::string pyne::nucname::openmc(std::string nuc) {
  return openmc(id(nuc));
}

//
// OPENMC -> id
//
int pyne::nucname::openmc_to_id(const char * nuc) {
  return openmc_to_id(std::string(nuc));
}

int pyne::nucname::openmc_to_id(std::string nuc) {
  std::string nucname;
  name_zz_t zznames = get_name_zz();

  // first two characters
  std::string::iterator aaa_start;
  int zzz = 0;
  if (zznames.count(nuc.substr(0,2)) == 1) {
    aaa_start = nuc.begin() + 2;
    zzz = zznames[nuc.substr(0,2)];
  }
  // then try only the first
  else if (zznames.count(nuc.substr(0,1)) == 1) {
    aaa_start = nuc.begin() + 1;
    zzz = zznames[nuc.substr(0,1)];
  } else {
    throw NotANuclide(nuc, "Not in the OpenMC format");
  }

  // set aaa - stop on "-" if the character exists
  std::string::iterator aaa_end = std::find(nuc.begin(), nuc.end(), '_');
  int aaa = pyne::to_int(nuc.substr(aaa_start - nuc.begin(), aaa_end - aaa_start));

  // check for metastable state
  int m = 0;
  if (aaa_end != nuc.end()) {
    std::string::iterator m_start = aaa_end + 2; // move forward once to skip "_m" characters
    m = pyne::to_int(nuc.substr(m_start - nuc.begin(), nuc.end() - m_start));
  }

  // form integer id and return
  return (zzz * 10000000) + (aaa * 10000) + m;

}


/**********************/
/*** fluka functions ***/
/**********************/
std::string pyne::nucname::fluka(int nuc) {
  int x = id(nuc);
  if (zz_fluka.count(x) == 0) {
    throw NotANuclide(nuc, "fluka name could not be found");
  }
  return zz_fluka[x];
}


//
// FLUKA name -> id
//
int pyne::nucname::fluka_to_id(std::string name) {
  if (fluka_zz.count(name) == 0) {
    throw NotANuclide(-1, "No nuclide: fluka name could not be found");
  }
  return fluka_zz[name];
}

int pyne::nucname::fluka_to_id(char * name) {
  return fluka_to_id(std::string(name));
}


/*************************/
/*** serpent functions ***/
/*************************/
std::string pyne::nucname::serpent(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";

  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int zzz = nucid / 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add LL
  std::string llupper = pyne::to_upper(zz_name[zzz]);
  std::string lllower = pyne::to_lower(zz_name[zzz]);
  newnuc += llupper[0];
  for (int l = 1; l < lllower.size(); l++)
    newnuc += lllower[l];

  // Add required dash
  newnuc += "-";

  // Add A-number
  if (0 < aaassss)
    newnuc += pyne::to_str(aaa);
  else if (0 == aaassss)
    newnuc += "nat";

  // Add meta-stable flag
  if (0 < ssss)
    newnuc += "m";

  return newnuc;
}


std::string pyne::nucname::serpent(const char * nuc) {
  std::string newnuc (nuc);
  return serpent(newnuc);
}


std::string pyne::nucname::serpent(std::string nuc) {
  return serpent(id(nuc));
}

//
// Serpent -> id
//
//int pyne::nucname::serpent_to_id(int nuc)
//{
// Should be ZAID
//}


int pyne::nucname::serpent_to_id(const char * nuc) {
  return serpent_to_id(std::string(nuc));
}


int pyne::nucname::serpent_to_id(std::string nuc) {
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int nucid;
  std::string elem_name;

  // Get the string into a regular form
  std::string nucstr = pyne::to_upper(nuc);
  nucstr = pyne::remove_substring(nucstr, "-");
  int nuclen = nucstr.length();

  // Nuclide is probably in name form, or some variation therein
  std::string anum_str = pyne::remove_characters(nucstr, pyne::alphabet);

  // natural element form, a la 'U' -> 920000000
  if (anum_str.empty() || pyne::contains_substring(nucstr, "NAT")) {
    elem_name = pyne::capitalize(pyne::remove_substring(nucstr, "NAT"));
    if (0 < name_zz.count(elem_name))
      return 10000000 * name_zz[elem_name];
  }
  int anum = pyne::to_int(anum_str);

  // Figure out if we are meta-stable or not
  std::string end_char = pyne::last_char(nucstr);
  if (end_char == "M")
    nucid = (10000 * anum) + 1;
  else if (pyne::contains_substring(pyne::digits, end_char))
    nucid = (10000 * anum);
  else
    throw NotANuclide(nucstr, nucid);

  // Add the Z-number
  elem_name = pyne::remove_characters(nucstr.substr(0, nuclen-1), pyne::digits);
  elem_name = pyne::capitalize(elem_name);
  if (0 < name_zz.count(elem_name))
    nucid = (10000000 * name_zz[elem_name]) + nucid;
  else
    throw NotANuclide(nucstr, nucid);
  return nucid;
}


/**********************/
/*** nist functions ***/
/**********************/
std::string pyne::nucname::nist(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add A-number
  if (0 < aaassss)
    newnuc += pyne::to_str(aaa);

  // Add name
  std::string name_upper = pyne::to_upper(zz_name[zzz]);
  std::string name_lower = pyne::to_lower(zz_name[zzz]);
  newnuc += name_upper[0];
  for (int l = 1; l < name_lower.size(); l++)
    newnuc += name_lower[l];

  // Add meta-stable flag
  // No metastable flag for NIST,
  // but could add star, by uncommenting below
  //if (0 < mod_10)
  //  newnuc += "*";

  return newnuc;
}


std::string pyne::nucname::nist(const char * nuc) {
  std::string newnuc (nuc);
  return nist(newnuc);
}


std::string pyne::nucname::nist(std::string nuc) {
  return nist(id(nuc));
}


//
// NIST -> id
//
//int pyne::nucname::nist_to_id(int nuc)
//{
// NON-EXISTANT
//};

int pyne::nucname::nist_to_id(const char * nuc) {
  return nist_to_id(std::string(nuc));
}

int pyne::nucname::nist_to_id(std::string nuc) {
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int nucid;
  nuc = pyne::to_upper(nuc);
  std::string elem_name;
  int nuclen = nuc.length();

  // Nuclide is probably in name form, or some variation therein
  std::string anum_str = pyne::remove_characters(nuc, pyne::alphabet);

  // natural element form, a la 'U' -> 920000000
  if (anum_str.empty()) {
    elem_name = pyne::capitalize(nuc);
    if (0 < name_zz.count(elem_name))
      return 10000000 * name_zz[elem_name];
  }
  nucid = pyne::to_int(anum_str) * 10000;

  // Add the Z-number
  elem_name = pyne::remove_characters(nuc, pyne::digits);
  elem_name = pyne::capitalize(elem_name);
  if (0 < name_zz.count(elem_name))
    nucid = (10000000 * name_zz[elem_name]) + nucid;
  else
    throw NotANuclide(nuc, nucid);
  return nucid;
}


/************************/
/*** cinder functions ***/
/************************/
int pyne::nucname::cinder(int nuc) {
  // cinder nuclides of form aaazzzm
  int nucid = id(nuc);
  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;
  if (10 <= ssss)
    ssss = 9;
  return (aaa*10000) + (zzz*10) + ssss;
}



int pyne::nucname::cinder(const char * nuc) {
  std::string newnuc (nuc);
  return cinder(newnuc);
}



int pyne::nucname::cinder(std::string nuc) {
  return cinder(id(nuc));
}

//
// Cinder -> Id
//
int pyne::nucname::cinder_to_id(int nuc) {
  int ssss = nuc % 10;
  int aaazzz = nuc / 10;
  int zzz = aaazzz % 1000;
  int aaa = aaazzz / 1000;
  return (zzz * 10000000) + (aaa * 10000) + ssss;
}


int pyne::nucname::cinder_to_id(const char * nuc) {
  return cinder_to_id(std::string(nuc));
}


int pyne::nucname::cinder_to_id(std::string nuc) {
  return cinder_to_id(pyne::to_int(nuc));
}




/**********************/
/*** ALARA functions ***/
/**********************/
std::string pyne::nucname::alara(int nuc) {
  int nucid = id(nuc);
  std::string newnuc = "";
  std::string ll = "";

  int zzz = nucid / 10000000;
  int ssss = nucid % 10000;
  int aaassss = nucid % 10000000;
  int aaa = aaassss / 10000;

  // Make sure the LL value is correct
  if (0 == zz_name.count(zzz))
    throw NotANuclide(nuc, nucid);

  // Add LL, in lower case
  ll += zz_name[zzz];

  for(int i = 0; ll[i] != '\0'; i++)
    ll[i] = tolower(ll[i]);
  newnuc += ll;

  // Add A-number
  if (0 < aaassss){
    newnuc += ":";
    newnuc += pyne::to_str(aaa);
  }

  // Note, ALARA input format does not use metastable flag
  return newnuc;
}


std::string pyne::nucname::alara(const char * nuc) {
  std::string newnuc (nuc);
  return alara(newnuc);
}


std::string pyne::nucname::alara(std::string nuc) {
  return alara(id(nuc));
}


//
// Cinder -> Id
//
//int pyne::nucname::alara_to_id(int nuc)
//{
// Not Possible
//}


int pyne::nucname::alara_to_id(const char * nuc) {
  return alara_to_id(std::string(nuc));
}


int pyne::nucname::alara_to_id(std::string nuc) {
  if (nuc.empty())
    throw NotANuclide(nuc, "<empty>");
  int nucid;
  nuc = pyne::to_upper(pyne::remove_characters(nuc, ":"));
  std::string elem_name;
  int nuclen = nuc.length();

  // Nuclide is probably in name form, or some variation therein
  std::string anum_str = pyne::remove_characters(nuc, pyne::alphabet);

  // natural element form, a la 'U' -> 920000000
  if (anum_str.empty()) {
    elem_name = pyne::capitalize(nuc);
    if (0 < name_zz.count(elem_name))
      return 10000000 * name_zz[elem_name];
  }
  nucid = pyne::to_int(anum_str) * 10000;

  // Add the Z-number
  elem_name = pyne::remove_characters(nuc, pyne::digits);
  elem_name = pyne::capitalize(elem_name);
  if (0 < name_zz.count(elem_name))
    nucid = (10000000 * name_zz[elem_name]) + nucid;
  else
    throw NotANuclide(nuc, nucid);
  return nucid;
}




/***********************/
/***  SZA functions  ***/
/***********************/
int pyne::nucname::sza(int nuc) {
  int nucid = id(nuc);
  int zzzaaa = nucid / 10000;
  int sss = nucid % 10000;
  return sss * 1000000 + zzzaaa;
}


int pyne::nucname::sza(const char * nuc) {
  std::string newnuc (nuc);
  return sza(newnuc);
}


int pyne::nucname::sza(std::string nuc) {
  return sza(id(nuc));
}


int pyne::nucname::sza_to_id(int nuc) {
  int sss = nuc / 1000000;
  int zzzaaa = nuc % 1000000;
  if (5 < sss){
  // Unphysical metastable state warning
   warning("You have indicated a metastable state of " + pyne::to_str(sss) + ". Metastable state above 5, possibly unphysical. ");
  }
  return zzzaaa * 10000 + sss;
}


int pyne::nucname::sza_to_id(const char * nuc) {
  std::string newnuc (nuc);
  return sza_to_id(newnuc);
}


int pyne::nucname::sza_to_id(std::string nuc) {
  return sza_to_id(pyne::to_int(nuc));
}


void pyne::nucname::_load_state_map(){
    for (int i = 0; i < TOTAL_STATE_MAPS; ++i) {
       state_id_map[map_nuc_ids[i]] = map_metastable[i];
    }
}

int pyne::nucname::state_id_to_id(int state) {
  int zzzaaa = (state / 10000) * 10000;
  int state_number = state % 10000;
  if (state_number == 0) return state;
  std::map<int, int>::iterator nuc_iter, nuc_end;

  nuc_iter = state_id_map.find(state);
  nuc_end = state_id_map.end();
  if (nuc_iter != nuc_end){
    int m = (*nuc_iter).second;
    return zzzaaa + m;
  }

  if (state_id_map.empty())  {
    _load_state_map();
    return state_id_to_id(state);
  }
  return -1;
}


int pyne::nucname::id_to_state_id(int nuc_id) {
  int zzzaaa = (nuc_id / 10000) * 10000;
  int state = nuc_id % 10000;
  if (state == 0) return nuc_id;
  std::map<int, int>::iterator nuc_iter, nuc_end, it;

  nuc_iter = state_id_map.lower_bound(zzzaaa);
  nuc_end = state_id_map.upper_bound(zzzaaa + 9999);
  for (it = nuc_iter; it!= nuc_end; ++it){
    if (state == it->second) {
      return it->first;
    }
  }

  if (state_id_map.empty())  {
    _load_state_map();
    return id_to_state_id(nuc_id);
  }
  return -1;
}


/************************/
/*** ENSDF functions ***/
/************************/
//
// ENSDF  -> Id
//

int pyne::nucname::ensdf_to_id(const char * nuc) {
  return ensdf_to_id(std::string(nuc));
}

int pyne::nucname::ensdf_to_id(std::string nuc) {
  if (nuc.size() < 4) {
    return nucname::id(nuc);
  } else if (std::isdigit(nuc[3])) {
    int aaa = to_int(nuc.substr(0, 3));
    int zzz;
    std::string xx_str = nuc.substr(3,2);
    zzz = to_int(xx_str) + 100;
    int nid = 10000 * aaa + 10000000 * zzz;
    return nid;
  } else {
    return nucname::id(nuc);
  }

}

//
// end of /home/mouginot/work/app/pyne/src/nucname.cpp
//


//
// start of /home/mouginot/work/app/pyne/src/jsoncpp.cpp
//
/// Json-cpp amalgated source (http://jsoncpp.sourceforge.net/).
/// It is intented to be used with #include <json.h>

// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: LICENSE
// //////////////////////////////////////////////////////////////////////

/*
The JsonCpp library's source code, including accompanying documentation,
tests and demonstration applications, are licensed under the following
conditions...

The author (Baptiste Lepilleur) explicitly disclaims copyright in all
jurisdictions which recognize such a disclaimer. In such jurisdictions,
this software is released into the Public Domain.

In jurisdictions which do not recognize Public Domain property (e.g. Germany as of
2010), this software is Copyright (c) 2007-2010 by Baptiste Lepilleur, and is
released under the terms of the MIT License (see below).

In jurisdictions which recognize Public Domain property, the user of this
software may choose to accept it either as 1) Public Domain, 2) under the
conditions of the MIT License (see below), or 3) under the terms of dual
Public Domain/MIT License conditions described here, as they choose.

The MIT License is about as close to Public Domain as a license can get, and is
described in clear, concise terms at:

   http://en.wikipedia.org/wiki/MIT_License

The full text of the MIT License follows:

========================================================================
Copyright (c) 2007-2010 Baptiste Lepilleur

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
========================================================================
(END LICENSE TEXT)

The MIT license is compatible with both the GPL and commercial
software, affording one all of the rights of Public Domain with the
minor nuisance of being required to keep the above copyright notice
and license text in the source code. Note also that by accepting the
Public Domain "license" you can re-license your copy using whatever
license you like.

*/

// //////////////////////////////////////////////////////////////////////
// End of content of file: LICENSE
// //////////////////////////////////////////////////////////////////////





#ifdef PYNE_IS_AMALGAMATED
  #if !defined(JSON_IS_AMALGAMATION)
    #define JSON_IS_AMALGAMATION
  #endif
#else
  #include "json.h"
#endif


// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_tool.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef LIB_JSONCPP_JSON_TOOL_H_INCLUDED
# define LIB_JSONCPP_JSON_TOOL_H_INCLUDED

/* This header provides common string manipulation support, such as UTF-8,
 * portable conversion from/to string...
 *
 * It is an internal header that must not be exposed.
 */

namespace Json {

/// Converts a unicode code-point to UTF-8.
static inline std::string
codePointToUTF8(unsigned int cp) {
   std::string result;

   // based on description from http://en.wikipedia.org/wiki/UTF-8

   if (cp <= 0x7f) {
      result.resize(1);
      result[0] = static_cast<char>(cp);
   } else if (cp <= 0x7FF) {
      result.resize(2);
      result[1] = static_cast<char>(0x80 | (0x3f & cp));
      result[0] = static_cast<char>(0xC0 | (0x1f & (cp >> 6)));
   } else if (cp <= 0xFFFF) {
      result.resize(3);
      result[2] = static_cast<char>(0x80 | (0x3f & cp));
      result[1] = 0x80 | static_cast<char>((0x3f & (cp >> 6)));
      result[0] = 0xE0 | static_cast<char>((0xf & (cp >> 12)));
   } else if (cp <= 0x10FFFF) {
      result.resize(4);
      result[3] = static_cast<char>(0x80 | (0x3f & cp));
      result[2] = static_cast<char>(0x80 | (0x3f & (cp >> 6)));
      result[1] = static_cast<char>(0x80 | (0x3f & (cp >> 12)));
      result[0] = static_cast<char>(0xF0 | (0x7 & (cp >> 18)));
   }

   return result;
}


/// Returns true if ch is a control character (in range [0,32[).
static inline bool
isControlCharacter(char ch) {
   return ch > 0 && ch <= 0x1F;
}


enum {
   /// Constant that specify the size of the buffer that must be passed to uintToString.
   uintToStringBufferSize = 3*sizeof(LargestUInt)+1
};

// Defines a char buffer for use with uintToString().
typedef char UIntToStringBuffer[uintToStringBufferSize];


/** Converts an unsigned integer to string.
 * @param value Unsigned interger to convert to string
 * @param current Input/Output string buffer.
 *        Must have at least uintToStringBufferSize chars free.
 */
static inline void
uintToString( LargestUInt value,
              char *&current ) {
   *--current = 0;
   do {
      *--current = char(value % 10) + '0';
      value /= 10;
   }
   while ( value != 0 );
}

} // namespace Json {

#endif // LIB_JSONCPP_JSON_TOOL_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_tool.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_reader.cpp
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#if !defined(JSON_IS_AMALGAMATION)
# include <json/reader.h>
# include <json/value.h>
# include "json_tool.h"
#endif // if !defined(JSON_IS_AMALGAMATION)
#include <utility>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <iostream>
#include <stdexcept>

#if _MSC_VER >= 1400 // VC++ 8.0
#pragma warning( disable : 4996 )   // disable warning about strdup being deprecated.
#endif

namespace Json {

// Implementation of class Features
// ////////////////////////////////

Features::Features()
   : allowComments_( true )
   , strictRoot_( false ) {
}


Features
Features::all() {
   return Features();
}


Features
Features::strictMode() {
   Features features;
   features.allowComments_ = false;
   features.strictRoot_ = true;
   return features;
}

// Implementation of class Reader
// ////////////////////////////////


static inline bool
in( Reader::Char c, Reader::Char c1, Reader::Char c2, Reader::Char c3, Reader::Char c4 ) {
   return c == c1  ||  c == c2  ||  c == c3  ||  c == c4;
}

static inline bool
in( Reader::Char c, Reader::Char c1, Reader::Char c2, Reader::Char c3, Reader::Char c4, Reader::Char c5 ) {
   return c == c1  ||  c == c2  ||  c == c3  ||  c == c4  ||  c == c5;
}


static bool
containsNewLine( Reader::Location begin,
                 Reader::Location end ) {
   for ( ;begin < end; ++begin )
      if ( *begin == '\n'  ||  *begin == '\r' )
         return true;
   return false;
}


// Class Reader
// //////////////////////////////////////////////////////////////////

Reader::Reader()
   : features_( Features::all() ) {
}


Reader::Reader( const Features &features )
   : features_( features ) {
}


bool
Reader::parse( const std::string &document,
               Value &root,
               bool collectComments ) {
   document_ = document;
   const char *begin = document_.c_str();
   const char *end = begin + document_.length();
   return parse( begin, end, root, collectComments );
}


bool
Reader::parse( std::istream& sin,
               Value &root,
               bool collectComments ) {
   //std::istream_iterator<char> begin(sin);
   //std::istream_iterator<char> end;
   // Those would allow streamed input from a file, if parse() were a
   // template function.

   // Since std::string is reference-counted, this at least does not
   // create an extra copy.
   std::string doc;
   std::getline(sin, doc, (char)EOF);
   return parse( doc, root, collectComments );
}

bool
Reader::parse( const char *beginDoc, const char *endDoc,
               Value &root,
               bool collectComments ) {
   if ( !features_.allowComments_ ) {
      collectComments = false;
   }

   begin_ = beginDoc;
   end_ = endDoc;
   collectComments_ = collectComments;
   current_ = begin_;
   lastValueEnd_ = 0;
   lastValue_ = 0;
   commentsBefore_ = "";
   errors_.clear();
   while ( !nodes_.empty() )
      nodes_.pop();
   nodes_.push( &root );

   bool successful = readValue();
   Token token;
   skipCommentTokens( token );
   if ( collectComments_  &&  !commentsBefore_.empty() )
      root.setComment( commentsBefore_, commentAfter );
   if ( features_.strictRoot_ ) {
      if ( !root.isArray()  &&  !root.isObject() ) {
         // Set error location to start of doc, ideally should be first token found in doc
         token.type_ = tokenError;
         token.start_ = beginDoc;
         token.end_ = endDoc;
         addError( "A valid JSON document must be either an array or an object value.",
                   token );
         return false;
      }
   }
   return successful;
}


bool
Reader::readValue() {
   Token token;
   skipCommentTokens( token );
   bool successful = true;

   if ( collectComments_  &&  !commentsBefore_.empty() ) {
      currentValue().setComment( commentsBefore_, commentBefore );
      commentsBefore_ = "";
   }


   switch ( token.type_ ) {
   case tokenObjectBegin:
      successful = readObject( token );
      break;
   case tokenArrayBegin:
      successful = readArray( token );
      break;
   case tokenNumber:
      successful = decodeNumber( token );
      break;
   case tokenString:
      successful = decodeString( token );
      break;
   case tokenTrue:
      currentValue() = true;
      break;
   case tokenFalse:
      currentValue() = false;
      break;
   case tokenNull:
      currentValue() = Value();
      break;
   default:
      return addError( "Syntax error: value, object or array expected.", token );
   }

   if ( collectComments_ ) {
      lastValueEnd_ = current_;
      lastValue_ = &currentValue();
   }

   return successful;
}


void
Reader::skipCommentTokens( Token &token ) {
   if ( features_.allowComments_ ) {
      do {
         readToken( token );
      }
      while ( token.type_ == tokenComment );
   } else {
      readToken( token );
   }
}


bool
Reader::expectToken( TokenType type, Token &token, const char *message ) {
   readToken( token );
   if ( token.type_ != type )
      return addError( message, token );
   return true;
}


bool
Reader::readToken( Token &token ) {
   skipSpaces();
   token.start_ = current_;
   Char c = getNextChar();
   bool ok = true;
   switch ( c ) {
   case '{':
      token.type_ = tokenObjectBegin;
      break;
   case '}':
      token.type_ = tokenObjectEnd;
      break;
   case '[':
      token.type_ = tokenArrayBegin;
      break;
   case ']':
      token.type_ = tokenArrayEnd;
      break;
   case '"':
      token.type_ = tokenString;
      ok = readString();
      break;
   case '/':
      token.type_ = tokenComment;
      ok = readComment();
      break;
   case '0':
   case '1':
   case '2':
   case '3':
   case '4':
   case '5':
   case '6':
   case '7':
   case '8':
   case '9':
   case '-':
      token.type_ = tokenNumber;
      readNumber();
      break;
   case 't':
      token.type_ = tokenTrue;
      ok = match( "rue", 3 );
      break;
   case 'f':
      token.type_ = tokenFalse;
      ok = match( "alse", 4 );
      break;
   case 'n':
      token.type_ = tokenNull;
      ok = match( "ull", 3 );
      break;
   case ',':
      token.type_ = tokenArraySeparator;
      break;
   case ':':
      token.type_ = tokenMemberSeparator;
      break;
   case 0:
      token.type_ = tokenEndOfStream;
      break;
   default:
      ok = false;
      break;
   }
   if ( !ok )
      token.type_ = tokenError;
   token.end_ = current_;
   return true;
}


void
Reader::skipSpaces() {
   while ( current_ != end_ ) {
      Char c = *current_;
      if ( c == ' '  ||  c == '\t'  ||  c == '\r'  ||  c == '\n' )
         ++current_;
      else
         break;
   }
}


bool
Reader::match( Location pattern,
               int patternLength ) {
   if ( end_ - current_ < patternLength )
      return false;
   int index = patternLength;
   while ( index-- )
      if ( current_[index] != pattern[index] )
         return false;
   current_ += patternLength;
   return true;
}


bool
Reader::readComment() {
   Location commentBegin = current_ - 1;
   Char c = getNextChar();
   bool successful = false;
   if ( c == '*' )
      successful = readCStyleComment();
   else if ( c == '/' )
      successful = readCppStyleComment();
   if ( !successful )
      return false;

   if ( collectComments_ ) {
      CommentPlacement placement = commentBefore;
      if ( lastValueEnd_  &&  !containsNewLine( lastValueEnd_, commentBegin ) )
      {
         if ( c != '*'  ||  !containsNewLine( commentBegin, current_ ) )
            placement = commentAfterOnSameLine;
      }

      addComment( commentBegin, current_, placement );
   }
   return true;
}


void
Reader::addComment( Location begin,
                    Location end,
                    CommentPlacement placement ) {
   assert( collectComments_ );
   if ( placement == commentAfterOnSameLine )
   {
      assert( lastValue_ != 0 );
      lastValue_->setComment( std::string( begin, end ), placement );
   } else {
      if ( !commentsBefore_.empty() )
         commentsBefore_ += "\n";
      commentsBefore_ += std::string( begin, end );
   }
}


bool
Reader::readCStyleComment() {
   while ( current_ != end_ ) {
      Char c = getNextChar();
      if ( c == '*'  &&  *current_ == '/' )
         break;
   }
   return getNextChar() == '/';
}


bool
Reader::readCppStyleComment() {
   while ( current_ != end_ ) {
      Char c = getNextChar();
      if (  c == '\r'  ||  c == '\n' )
         break;
   }
   return true;
}


void
Reader::readNumber() {
   while ( current_ != end_ ) {
      if ( !(*current_ >= '0'  &&  *current_ <= '9')  &&
           !in( *current_, '.', 'e', 'E', '+', '-' ) )
         break;
      ++current_;
   }
}

bool
Reader::readString() {
   Char c = 0;
   while ( current_ != end_ ) {
      c = getNextChar();
      if ( c == '\\' )
         getNextChar();
      else if ( c == '"' )
         break;
   }
   return c == '"';
}


bool
Reader::readObject( Token &/*tokenStart*/ ) {
   Token tokenName;
   std::string name;
   currentValue() = Value( objectValue );
   while ( readToken( tokenName ) ) {
      bool initialTokenOk = true;
      while ( tokenName.type_ == tokenComment  &&  initialTokenOk )
         initialTokenOk = readToken( tokenName );
      if  ( !initialTokenOk )
         break;
      if ( tokenName.type_ == tokenObjectEnd  &&  name.empty() )  // empty object
         return true;
      if ( tokenName.type_ != tokenString )
         break;

      name = "";
      if ( !decodeString( tokenName, name ) )
         return recoverFromError( tokenObjectEnd );

      Token colon;
      if ( !readToken( colon ) ||  colon.type_ != tokenMemberSeparator ) {
         return addErrorAndRecover( "Missing ':' after object member name",
                                    colon,
                                    tokenObjectEnd );
      }
      Value &value = currentValue()[ name ];
      nodes_.push( &value );
      bool ok = readValue();
      nodes_.pop();
      if ( !ok ) // error already set
         return recoverFromError( tokenObjectEnd );

      Token comma;
      if ( !readToken( comma )
            ||  ( comma.type_ != tokenObjectEnd  &&
                  comma.type_ != tokenArraySeparator &&
                  comma.type_ != tokenComment ) ) {
         return addErrorAndRecover( "Missing ',' or '}' in object declaration",
                                    comma,
                                    tokenObjectEnd );
      }
      bool finalizeTokenOk = true;
      while ( comma.type_ == tokenComment &&
              finalizeTokenOk )
         finalizeTokenOk = readToken( comma );
      if ( comma.type_ == tokenObjectEnd )
         return true;
   }
   return addErrorAndRecover( "Missing '}' or object member name",
                              tokenName,
                              tokenObjectEnd );
}


bool
Reader::readArray( Token &/*tokenStart*/ ) {
   currentValue() = Value( arrayValue );
   skipSpaces();
   if ( *current_ == ']' ) { // empty array
      Token endArray;
      readToken( endArray );
      return true;
   }
   int index = 0;
   for (;;) {
      Value &value = currentValue()[ index++ ];
      nodes_.push( &value );
      bool ok = readValue();
      nodes_.pop();
      if ( !ok ) // error already set
         return recoverFromError( tokenArrayEnd );

      Token token;
      // Accept Comment after last item in the array.
      ok = readToken( token );
      while ( token.type_ == tokenComment  &&  ok ) {
         ok = readToken( token );
      }
      bool badTokenType = ( token.type_ != tokenArraySeparator  &&
                            token.type_ != tokenArrayEnd );
      if ( !ok  ||  badTokenType ) {
         return addErrorAndRecover( "Missing ',' or ']' in array declaration",
                                    token,
                                    tokenArrayEnd );
      }
      if ( token.type_ == tokenArrayEnd )
         break;
   }
   return true;
}


bool
Reader::decodeNumber( Token &token ) {
   bool isDouble = false;
   for ( Location inspect = token.start_; inspect != token.end_; ++inspect )
   {
      isDouble = isDouble
                 ||  in( *inspect, '.', 'e', 'E', '+' )
                 ||  ( *inspect == '-'  &&  inspect != token.start_ );
   }
   if ( isDouble )
      return decodeDouble( token );
   // Attempts to parse the number as an integer. If the number is
   // larger than the maximum supported value of an integer then
   // we decode the number as a double.
   Location current = token.start_;
   bool isNegative = *current == '-';
   if ( isNegative )
      ++current;
   Value::LargestUInt maxIntegerValue = isNegative ? Value::LargestUInt(-Value::minLargestInt)
                                                   : Value::maxLargestUInt;
   Value::LargestUInt threshold = maxIntegerValue / 10;
   Value::UInt lastDigitThreshold = Value::UInt( maxIntegerValue % 10 );
   assert( lastDigitThreshold >=0  &&  lastDigitThreshold <= 9 );
   Value::LargestUInt value = 0;
   while ( current < token.end_ ) {
      Char c = *current++;
      if ( c < '0'  ||  c > '9' )
         return addError( "'" + std::string( token.start_, token.end_ ) + "' is not a number.", token );
      Value::UInt digit(c - '0');
      if ( value >= threshold ) {
         // If the current digit is not the last one, or if it is
         // greater than the last digit of the maximum integer value,
         // the parse the number as a double.
         if ( current != token.end_  ||  digit > lastDigitThreshold ) {
            return decodeDouble( token );
         }
      }
      value = value * 10 + digit;
   }
   if ( isNegative )
      currentValue() = -Value::LargestInt( value );
   else if ( value <= Value::LargestUInt(Value::maxInt) )
      currentValue() = Value::LargestInt( value );
   else
      currentValue() = value;
   return true;
}


bool
Reader::decodeDouble( Token &token ) {
   double value = 0;
   const int bufferSize = 32;
   int count;
   int length = int(token.end_ - token.start_);
   if ( length <= bufferSize ) {
      Char buffer[bufferSize+1];
      memcpy( buffer, token.start_, length );
      buffer[length] = 0;
      count = sscanf( buffer, "%lf", &value );
   } else {
      std::string buffer( token.start_, token.end_ );
      count = sscanf( buffer.c_str(), "%lf", &value );
   }

   if ( count != 1 )
      return addError( "'" + std::string( token.start_, token.end_ ) + "' is not a number.", token );
   currentValue() = value;
   return true;
}


bool
Reader::decodeString( Token &token ) {
   std::string decoded;
   if ( !decodeString( token, decoded ) )
      return false;
   currentValue() = decoded;
   return true;
}


bool
Reader::decodeString( Token &token, std::string &decoded ) {
   decoded.reserve( token.end_ - token.start_ - 2 );
   Location current = token.start_ + 1; // skip '"'
   Location end = token.end_ - 1;      // do not include '"'
   while ( current != end ) {
      Char c = *current++;
      if ( c == '"' )
         break;
      else if ( c == '\\' ) {
         if ( current == end )
            return addError( "Empty escape sequence in string", token, current );
         Char escape = *current++;
         switch ( escape ) {
         case '"': decoded += '"'; break;
         case '/': decoded += '/'; break;
         case '\\': decoded += '\\'; break;
         case 'b': decoded += '\b'; break;
         case 'f': decoded += '\f'; break;
         case 'n': decoded += '\n'; break;
         case 'r': decoded += '\r'; break;
         case 't': decoded += '\t'; break;
         case 'u': {
               unsigned int unicode;
               if ( !decodeUnicodeCodePoint( token, current, end, unicode ) )
                  return false;
               decoded += codePointToUTF8(unicode);
            }
            break;
         default:
            return addError( "Bad escape sequence in string", token, current );
         }
      }
      else {
         decoded += c;
      }
   }
   return true;
}

bool
Reader::decodeUnicodeCodePoint( Token &token,
                                     Location &current,
                                     Location end,
                                     unsigned int &unicode ) {

   if ( !decodeUnicodeEscapeSequence( token, current, end, unicode ) )
      return false;
   if (unicode >= 0xD800 && unicode <= 0xDBFF) {
      // surrogate pairs
      if (end - current < 6)
         return addError( "additional six characters expected to parse unicode surrogate pair.", token, current );
      unsigned int surrogatePair;
      if (*(current++) == '\\' && *(current++)== 'u') {
         if (decodeUnicodeEscapeSequence( token, current, end, surrogatePair )) {
            unicode = 0x10000 + ((unicode & 0x3FF) << 10) + (surrogatePair & 0x3FF);
         } else
            return false;
      } else
         return addError( "expecting another \\u token to begin the second half of a unicode surrogate pair", token, current );
   }
   return true;
}

bool
Reader::decodeUnicodeEscapeSequence( Token &token,
                                     Location &current,
                                     Location end,
                                     unsigned int &unicode ) {
   if ( end - current < 4 )
      return addError( "Bad unicode escape sequence in string: four digits expected.", token, current );
   unicode = 0;
   for ( int index =0; index < 4; ++index ) {
      Char c = *current++;
      unicode *= 16;
      if ( c >= '0'  &&  c <= '9' )
         unicode += c - '0';
      else if ( c >= 'a'  &&  c <= 'f' )
         unicode += c - 'a' + 10;
      else if ( c >= 'A'  &&  c <= 'F' )
         unicode += c - 'A' + 10;
      else
         return addError( "Bad unicode escape sequence in string: hexadecimal digit expected.", token, current );
   }
   return true;
}


bool
Reader::addError( const std::string &message,
                  Token &token,
                  Location extra ) {
   ErrorInfo info;
   info.token_ = token;
   info.message_ = message;
   info.extra_ = extra;
   errors_.push_back( info );
   return false;
}


bool
Reader::recoverFromError( TokenType skipUntilToken ) {
   int errorCount = int(errors_.size());
   Token skip;
   for (;;) {
      if ( !readToken(skip) )
         errors_.resize( errorCount ); // discard errors caused by recovery
      if ( skip.type_ == skipUntilToken  ||  skip.type_ == tokenEndOfStream )
         break;
   }
   errors_.resize( errorCount );
   return false;
}


bool
Reader::addErrorAndRecover( const std::string &message,
                            Token &token,
                            TokenType skipUntilToken ) {
   addError( message, token );
   return recoverFromError( skipUntilToken );
}


Value &
Reader::currentValue() {
   return *(nodes_.top());
}


Reader::Char
Reader::getNextChar() {
   if ( current_ == end_ )
      return 0;
   return *current_++;
}


void
Reader::getLocationLineAndColumn( Location location,
                                  int &line,
                                  int &column ) const {
   Location current = begin_;
   Location lastLineStart = current;
   line = 0;
   while ( current < location  &&  current != end_ ) {
      Char c = *current++;
      if ( c == '\r' ) {
         if ( *current == '\n' )
            ++current;
         lastLineStart = current;
         ++line;
      } else if ( c == '\n' ) {
         lastLineStart = current;
         ++line;
      }
   }
   // column & line start at 1
   column = int(location - lastLineStart) + 1;
   ++line;
}


std::string
Reader::getLocationLineAndColumn( Location location ) const {
   int line, column;
   getLocationLineAndColumn( location, line, column );
   char buffer[18+16+16+1];
   sprintf( buffer, "Line %d, Column %d", line, column );
   return buffer;
}


// Deprecated. Preserved for backward compatibility
std::string
Reader::getFormatedErrorMessages() const {
    return getFormattedErrorMessages();
}


std::string
Reader::getFormattedErrorMessages() const {
   std::string formattedMessage;
   for ( Errors::const_iterator itError = errors_.begin();
         itError != errors_.end();
         ++itError ) {
      const ErrorInfo &error = *itError;
      formattedMessage += "* " + getLocationLineAndColumn( error.token_.start_ ) + "\n";
      formattedMessage += "  " + error.message_ + "\n";
      if ( error.extra_ )
         formattedMessage += "See " + getLocationLineAndColumn( error.extra_ ) + " for detail.\n";
   }
   return formattedMessage;
}


std::istream& operator>>( std::istream &sin, Value &root ) {
    Json::Reader reader;
    bool ok = reader.parse(sin, root, true);
    //JSON_ASSERT( ok );
    if (!ok) throw std::runtime_error(reader.getFormattedErrorMessages());
    return sin;
}


} // namespace Json

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_reader.cpp
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_batchallocator.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef JSONCPP_BATCHALLOCATOR_H_INCLUDED
# define JSONCPP_BATCHALLOCATOR_H_INCLUDED

# include <stdlib.h>
# include <assert.h>

# ifndef JSONCPP_DOC_EXCLUDE_IMPLEMENTATION

namespace Json {

/* Fast memory allocator.
 *
 * This memory allocator allocates memory for a batch of object (specified by
 * the page size, the number of object in each page).
 *
 * It does not allow the destruction of a single object. All the allocated objects
 * can be destroyed at once. The memory can be either released or reused for future
 * allocation.
 *
 * The in-place new operator must be used to construct the object using the pointer
 * returned by allocate.
 */
template<typename AllocatedType
        ,const unsigned int objectPerAllocation>
class BatchAllocator {
public:
   typedef AllocatedType Type;

   BatchAllocator( unsigned int objectsPerPage = 255 )
      : freeHead_( 0 )
      , objectsPerPage_( objectsPerPage ) {
//      printf( "Size: %d => %s\n", sizeof(AllocatedType), typeid(AllocatedType).name() );
      assert( sizeof(AllocatedType) * objectPerAllocation >= sizeof(AllocatedType *) ); // We must be able to store a slist in the object free space.
      assert( objectsPerPage >= 16 );
      batches_ = allocateBatch( 0 );   // allocated a dummy page
      currentBatch_ = batches_;
   }

   ~BatchAllocator() {
      for ( BatchInfo *batch = batches_; batch;  ) {
         BatchInfo *nextBatch = batch->next_;
         free( batch );
         batch = nextBatch;
      }
   }

   /// allocate space for an array of objectPerAllocation object.
   /// @warning it is the responsability of the caller to call objects constructors.
   AllocatedType *allocate() {
      if ( freeHead_ ) { // returns node from free list.
         AllocatedType *object = freeHead_;
         freeHead_ = *(AllocatedType **)object;
         return object;
      }
      if ( currentBatch_->used_ == currentBatch_->end_ ) {
         currentBatch_ = currentBatch_->next_;
         while ( currentBatch_  &&  currentBatch_->used_ == currentBatch_->end_ )
            currentBatch_ = currentBatch_->next_;

         if ( !currentBatch_  ) { // no free batch found, allocate a new one
            currentBatch_ = allocateBatch( objectsPerPage_ );
            currentBatch_->next_ = batches_; // insert at the head of the list
            batches_ = currentBatch_;
         }
      }
      AllocatedType *allocated = currentBatch_->used_;
      currentBatch_->used_ += objectPerAllocation;
      return allocated;
   }

   /// Release the object.
   /// @warning it is the responsability of the caller to actually destruct the object.
   void release( AllocatedType *object ) {
      assert( object != 0 );
      *(AllocatedType **)object = freeHead_;
      freeHead_ = object;
   }

private:
   struct BatchInfo {
      BatchInfo *next_;
      AllocatedType *used_;
      AllocatedType *end_;
      AllocatedType buffer_[objectPerAllocation];
   };

   // disabled copy constructor and assignement operator.
   BatchAllocator( const BatchAllocator & );
   void operator =( const BatchAllocator &);

   static BatchInfo *allocateBatch( unsigned int objectsPerPage ) {
      const unsigned int mallocSize = sizeof(BatchInfo) - sizeof(AllocatedType)* objectPerAllocation
                                + sizeof(AllocatedType) * objectPerAllocation * objectsPerPage;
      BatchInfo *batch = static_cast<BatchInfo*>( malloc( mallocSize ) );
      batch->next_ = 0;
      batch->used_ = batch->buffer_;
      batch->end_ = batch->buffer_ + objectsPerPage;
      return batch;
   }

   BatchInfo *batches_;
   BatchInfo *currentBatch_;
   /// Head of a single linked list within the allocated space of freeed object
   AllocatedType *freeHead_;
   unsigned int objectsPerPage_;
};


} // namespace Json

# endif // ifndef JSONCPP_DOC_INCLUDE_IMPLEMENTATION

#endif // JSONCPP_BATCHALLOCATOR_H_INCLUDED


// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_batchallocator.h
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_valueiterator.inl
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

// included by json_value.cpp

namespace Json {

// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class ValueIteratorBase
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////

ValueIteratorBase::ValueIteratorBase()
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   : current_()
   , isNull_( true ) {
}
#else
   : isArray_( true )
   , isNull_( true ) {
   iterator_.array_ = ValueInternalArray::IteratorState();
}
#endif


#ifndef JSON_VALUE_USE_INTERNAL_MAP
ValueIteratorBase::ValueIteratorBase( const Value::ObjectValues::iterator &current )
   : current_( current )
   , isNull_( false ) {
}
#else
ValueIteratorBase::ValueIteratorBase( const ValueInternalArray::IteratorState &state )
   : isArray_( true ) {
   iterator_.array_ = state;
}


ValueIteratorBase::ValueIteratorBase( const ValueInternalMap::IteratorState &state )
   : isArray_( false ) {
   iterator_.map_ = state;
}
#endif

Value &
ValueIteratorBase::deref() const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   return current_->second;
#else
   if ( isArray_ )
      return ValueInternalArray::dereference( iterator_.array_ );
   return ValueInternalMap::value( iterator_.map_ );
#endif
}


void
ValueIteratorBase::increment() {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   ++current_;
#else
   if ( isArray_ )
      ValueInternalArray::increment( iterator_.array_ );
   ValueInternalMap::increment( iterator_.map_ );
#endif
}


void
ValueIteratorBase::decrement() {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   --current_;
#else
   if ( isArray_ )
      ValueInternalArray::decrement( iterator_.array_ );
   ValueInternalMap::decrement( iterator_.map_ );
#endif
}


ValueIteratorBase::difference_type
ValueIteratorBase::computeDistance( const SelfType &other ) const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
# ifdef JSON_USE_CPPTL_SMALLMAP
   return current_ - other.current_;
# else
   // Iterator for null value are initialized using the default
   // constructor, which initialize current_ to the default
   // std::map::iterator. As begin() and end() are two instance
   // of the default std::map::iterator, they can not be compared.
   // To allow this, we handle this comparison specifically.
   if ( isNull_  &&  other.isNull_ ) {
      return 0;
   }


   // Usage of std::distance is not portable (does not compile with Sun Studio 12 RogueWave STL,
   // which is the one used by default).
   // Using a portable hand-made version for non random iterator instead:
   //   return difference_type( std::distance( current_, other.current_ ) );
   difference_type myDistance = 0;
   for ( Value::ObjectValues::iterator it = current_; it != other.current_; ++it ) {
      ++myDistance;
   }
   return myDistance;
# endif
#else
   if ( isArray_ )
      return ValueInternalArray::distance( iterator_.array_, other.iterator_.array_ );
   return ValueInternalMap::distance( iterator_.map_, other.iterator_.map_ );
#endif
}


bool
ValueIteratorBase::isEqual( const SelfType &other ) const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   if ( isNull_ ) {
      return other.isNull_;
   }
   return current_ == other.current_;
#else
   if ( isArray_ )
      return ValueInternalArray::equals( iterator_.array_, other.iterator_.array_ );
   return ValueInternalMap::equals( iterator_.map_, other.iterator_.map_ );
#endif
}


void
ValueIteratorBase::copy( const SelfType &other ) {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   current_ = other.current_;
#else
   if ( isArray_ )
      iterator_.array_ = other.iterator_.array_;
   iterator_.map_ = other.iterator_.map_;
#endif
}


Value
ValueIteratorBase::key() const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   const Value::CZString czstring = (*current_).first;
   if ( czstring.c_str() ) {
      if ( czstring.isStaticString() )
         return Value( StaticString( czstring.c_str() ) );
      return Value( czstring.c_str() );
   }
   return Value( czstring.index() );
#else
   if ( isArray_ )
      return Value( ValueInternalArray::indexOf( iterator_.array_ ) );
   bool isStatic;
   const char *memberName = ValueInternalMap::key( iterator_.map_, isStatic );
   if ( isStatic )
      return Value( StaticString( memberName ) );
   return Value( memberName );
#endif
}


UInt
ValueIteratorBase::index() const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   const Value::CZString czstring = (*current_).first;
   if ( !czstring.c_str() )
      return czstring.index();
   return Value::UInt( -1 );
#else
   if ( isArray_ )
      return Value::UInt( ValueInternalArray::indexOf( iterator_.array_ ) );
   return Value::UInt( -1 );
#endif
}


const char *
ValueIteratorBase::memberName() const {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   const char *name = (*current_).first.c_str();
   return name ? name : "";
#else
   if ( !isArray_ )
      return ValueInternalMap::key( iterator_.map_ );
   return "";
#endif
}


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class ValueConstIterator
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////

ValueConstIterator::ValueConstIterator() {
}


#ifndef JSON_VALUE_USE_INTERNAL_MAP
ValueConstIterator::ValueConstIterator( const Value::ObjectValues::iterator &current )
   : ValueIteratorBase( current ) {
}
#else
ValueConstIterator::ValueConstIterator( const ValueInternalArray::IteratorState &state )
   : ValueIteratorBase( state ) {
}

ValueConstIterator::ValueConstIterator( const ValueInternalMap::IteratorState &state )
   : ValueIteratorBase( state ) {
}
#endif

ValueConstIterator &
ValueConstIterator::operator =( const ValueIteratorBase &other ) {
   copy( other );
   return *this;
}


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class ValueIterator
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////

ValueIterator::ValueIterator() {
}


#ifndef JSON_VALUE_USE_INTERNAL_MAP
ValueIterator::ValueIterator( const Value::ObjectValues::iterator &current )
   : ValueIteratorBase( current ) {
}
#else
ValueIterator::ValueIterator( const ValueInternalArray::IteratorState &state )
   : ValueIteratorBase( state ) {
}

ValueIterator::ValueIterator( const ValueInternalMap::IteratorState &state )
   : ValueIteratorBase( state ) {
}
#endif

ValueIterator::ValueIterator( const ValueConstIterator &other )
   : ValueIteratorBase( other ) {
}

ValueIterator::ValueIterator( const ValueIterator &other )
   : ValueIteratorBase( other ) {
}

ValueIterator &
ValueIterator::operator =( const SelfType &other ) {
   copy( other );
   return *this;
}

} // namespace Json

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_valueiterator.inl
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_value.cpp
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#if !defined(JSON_IS_AMALGAMATION)
# include <json/value.h>
# include <json/writer.h>
# ifndef JSON_USE_SIMPLE_INTERNAL_ALLOCATOR
#  include "json_batchallocator.h"
# endif // #ifndef JSON_USE_SIMPLE_INTERNAL_ALLOCATOR
#endif // if !defined(JSON_IS_AMALGAMATION)
#include <iostream>
#include <utility>
#include <stdexcept>
#include <cstring>
#include <cassert>
#ifdef JSON_USE_CPPTL
# include <cpptl/conststring.h>
#endif
#include <cstddef>    // size_t

#define JSON_ASSERT_UNREACHABLE assert( false )
#define JSON_ASSERT( condition ) assert( condition );  // @todo <= change this into an exception throw
#define JSON_FAIL_MESSAGE( message ) throw std::runtime_error( message );
#define JSON_ASSERT_MESSAGE( condition, message ) if (!( condition )) JSON_FAIL_MESSAGE( message )

namespace Json {

const Value Value::null;
const Int Value::minInt = Int( ~(UInt(-1)/2) );
const Int Value::maxInt = Int( UInt(-1)/2 );
const UInt Value::maxUInt = UInt(-1);
const Int64 Value::minInt64 = Int64( ~(UInt64(-1)/2) );
const Int64 Value::maxInt64 = Int64( UInt64(-1)/2 );
const UInt64 Value::maxUInt64 = UInt64(-1);
const LargestInt Value::minLargestInt = LargestInt( ~(LargestUInt(-1)/2) );
const LargestInt Value::maxLargestInt = LargestInt( LargestUInt(-1)/2 );
const LargestUInt Value::maxLargestUInt = LargestUInt(-1);


/// Unknown size marker
static const unsigned int unknown = (unsigned)-1;


/** Duplicates the specified string value.
 * @param value Pointer to the string to duplicate. Must be zero-terminated if
 *              length is "unknown".
 * @param length Length of the value. if equals to unknown, then it will be
 *               computed using strlen(value).
 * @return Pointer on the duplicate instance of string.
 */
static inline char *
duplicateStringValue( const char *value,
                      unsigned int length = unknown ) {
   if ( length == unknown )
      length = (unsigned int)strlen(value);
   char *newString = static_cast<char *>( malloc( length + 1 ) );
   JSON_ASSERT_MESSAGE( newString != 0, "Failed to allocate string value buffer" );
   memcpy( newString, value, length );
   newString[length] = 0;
   return newString;
}


/** Free the string duplicated by duplicateStringValue().
 */
static inline void
releaseStringValue( char *value ) {
   if ( value )
      free( value );
}

} // namespace Json


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// ValueInternals...
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
#if !defined(JSON_IS_AMALGAMATION)
# ifdef JSON_VALUE_USE_INTERNAL_MAP
#  include "json_internalarray.inl"
#  include "json_internalmap.inl"
# endif // JSON_VALUE_USE_INTERNAL_MAP

# include "json_valueiterator.inl"
#endif // if !defined(JSON_IS_AMALGAMATION)

namespace Json {

// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class Value::CommentInfo
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////


Value::CommentInfo::CommentInfo()
   : comment_( 0 ) {
}

Value::CommentInfo::~CommentInfo() {
   if ( comment_ )
      releaseStringValue( comment_ );
}


void
Value::CommentInfo::setComment( const char *text ) {
   if ( comment_ )
      releaseStringValue( comment_ );
   JSON_ASSERT( text != 0 );
   JSON_ASSERT_MESSAGE( text[0]=='\0' || text[0]=='/', "Comments must start with /");
   // It seems that /**/ style comments are acceptable as well.
   comment_ = duplicateStringValue( text );
}


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class Value::CZString
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
# ifndef JSON_VALUE_USE_INTERNAL_MAP

// Notes: index_ indicates if the string was allocated when
// a string is stored.

Value::CZString::CZString( ArrayIndex index )
   : cstr_( 0 )
   , index_( index ) {
}

Value::CZString::CZString( const char *cstr, DuplicationPolicy allocate )
   : cstr_( allocate == duplicate ? duplicateStringValue(cstr)
                                  : cstr )
   , index_( allocate ) {
}

Value::CZString::CZString( const CZString &other )
: cstr_( other.index_ != noDuplication &&  other.cstr_ != 0
                ?  duplicateStringValue( other.cstr_ )
                : other.cstr_ )
   , index_( other.cstr_ ? (other.index_ == noDuplication ? noDuplication : duplicate)
                         : other.index_ ) {
}

Value::CZString::~CZString() {
   if ( cstr_  &&  index_ == duplicate )
      releaseStringValue( const_cast<char *>( cstr_ ) );
}

void
Value::CZString::swap( CZString &other ) {
   std::swap( cstr_, other.cstr_ );
   std::swap( index_, other.index_ );
}

Value::CZString &
Value::CZString::operator =( const CZString &other ) {
   CZString temp( other );
   swap( temp );
   return *this;
}

bool
Value::CZString::operator<( const CZString &other ) const  {
   if ( cstr_ )
      return strcmp( cstr_, other.cstr_ ) < 0;
   return index_ < other.index_;
}

bool
Value::CZString::operator==( const CZString &other ) const  {
   if ( cstr_ )
      return strcmp( cstr_, other.cstr_ ) == 0;
   return index_ == other.index_;
}


ArrayIndex
Value::CZString::index() const {
   return index_;
}


const char *
Value::CZString::c_str() const {
   return cstr_;
}

bool
Value::CZString::isStaticString() const {
   return index_ == noDuplication;
}

#endif // ifndef JSON_VALUE_USE_INTERNAL_MAP


// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// class Value::Value
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////

/*! \internal Default constructor initialization must be equivalent to:
 * memset( this, 0, sizeof(Value) )
 * This optimization is used in ValueInternalMap fast allocator.
 */
Value::Value( ValueType type )
   : type_( type )
   , allocated_( 0 )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   switch ( type ) {
   case nullValue:
      break;
   case intValue:
   case uintValue:
      value_.int_ = 0;
      break;
   case realValue:
      value_.real_ = 0.0;
      break;
   case stringValue:
      value_.string_ = 0;
      break;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
   case objectValue:
      value_.map_ = new ObjectValues();
      break;
#else
   case arrayValue:
      value_.array_ = arrayAllocator()->newArray();
      break;
   case objectValue:
      value_.map_ = mapAllocator()->newMap();
      break;
#endif
   case booleanValue:
      value_.bool_ = false;
      break;
   default:
      JSON_ASSERT_UNREACHABLE;
   }
}


#if defined(JSON_HAS_INT64)
Value::Value( UInt value )
   : type_( uintValue )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.uint_ = value;
}

Value::Value( Int value )
   : type_( intValue )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.int_ = value;
}

#endif // if defined(JSON_HAS_INT64)


Value::Value( Int64 value )
   : type_( intValue )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.int_ = value;
}


Value::Value( UInt64 value )
   : type_( uintValue )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.uint_ = value;
}

Value::Value( double value )
   : type_( realValue )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.real_ = value;
}

Value::Value( const char *value )
   : type_( stringValue )
   , allocated_( true )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.string_ = duplicateStringValue( value );
}


Value::Value( const char *beginValue,
              const char *endValue )
   : type_( stringValue )
   , allocated_( true )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.string_ = duplicateStringValue( beginValue,
                                          (unsigned int)(endValue - beginValue) );
}


Value::Value( const std::string &value )
   : type_( stringValue )
   , allocated_( true )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.string_ = duplicateStringValue( value.c_str(),
                                          (unsigned int)value.length() );

}

Value::Value( const StaticString &value )
   : type_( stringValue )
   , allocated_( false )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.string_ = const_cast<char *>( value.c_str() );
}


# ifdef JSON_USE_CPPTL
Value::Value( const CppTL::ConstString &value )
   : type_( stringValue )
   , allocated_( true )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.string_ = duplicateStringValue( value, value.length() );
}
# endif

Value::Value( bool value )
   : type_( booleanValue )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   value_.bool_ = value;
}


Value::Value( const Value &other )
   : type_( other.type_ )
   , comments_( 0 )
# ifdef JSON_VALUE_USE_INTERNAL_MAP
   , itemIsUsed_( 0 )
#endif
{
   switch ( type_ ) {
   case nullValue:
   case intValue:
   case uintValue:
   case realValue:
   case booleanValue:
      value_ = other.value_;
      break;
   case stringValue:
      if ( other.value_.string_ ) {
         value_.string_ = duplicateStringValue( other.value_.string_ );
         allocated_ = true;
      }
      else
         value_.string_ = 0;
      break;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
   case objectValue:
      value_.map_ = new ObjectValues( *other.value_.map_ );
      break;
#else
   case arrayValue:
      value_.array_ = arrayAllocator()->newArrayCopy( *other.value_.array_ );
      break;
   case objectValue:
      value_.map_ = mapAllocator()->newMapCopy( *other.value_.map_ );
      break;
#endif
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   if ( other.comments_ ) {
      comments_ = new CommentInfo[numberOfCommentPlacement];
      for ( int comment =0; comment < numberOfCommentPlacement; ++comment ) {
         const CommentInfo &otherComment = other.comments_[comment];
         if ( otherComment.comment_ )
            comments_[comment].setComment( otherComment.comment_ );
      }
   }
}


Value::~Value() {
   switch ( type_ ) {
   case nullValue:
   case intValue:
   case uintValue:
   case realValue:
   case booleanValue:
      break;
   case stringValue:
      if ( allocated_ )
         releaseStringValue( value_.string_ );
      break;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
   case objectValue:
      delete value_.map_;
      break;
#else
   case arrayValue:
      arrayAllocator()->destructArray( value_.array_ );
      break;
   case objectValue:
      mapAllocator()->destructMap( value_.map_ );
      break;
#endif
   default:
      JSON_ASSERT_UNREACHABLE;
   }

   if ( comments_ )
      delete[] comments_;
}

Value &
Value::operator=( const Value &other ) {
   Value temp( other );
   swap( temp );
   return *this;
}

void
Value::swap( Value &other ) {
   ValueType temp = type_;
   type_ = other.type_;
   other.type_ = temp;
   std::swap( value_, other.value_ );
   int temp2 = allocated_;
   allocated_ = other.allocated_;
   other.allocated_ = temp2;
}

ValueType
Value::type() const {
   return type_;
}


int
Value::compare( const Value &other ) const {
   if ( *this < other )
      return -1;
   if ( *this > other )
      return 1;
   return 0;
}


bool
Value::operator <( const Value &other ) const {
   int typeDelta = type_ - other.type_;
   if ( typeDelta )
      return typeDelta < 0 ? true : false;
   switch ( type_ ) {
   case nullValue:
      return false;
   case intValue:
      return value_.int_ < other.value_.int_;
   case uintValue:
      return value_.uint_ < other.value_.uint_;
   case realValue:
      return value_.real_ < other.value_.real_;
   case booleanValue:
      return value_.bool_ < other.value_.bool_;
   case stringValue:
      return ( value_.string_ == 0  &&  other.value_.string_ )
             || ( other.value_.string_
                  &&  value_.string_
                  && strcmp( value_.string_, other.value_.string_ ) < 0 );
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
   case objectValue: {
         int delta = int( value_.map_->size() - other.value_.map_->size() );
         if ( delta )
            return delta < 0;
         return (*value_.map_) < (*other.value_.map_);
      }
#else
   case arrayValue:
      return value_.array_->compare( *(other.value_.array_) ) < 0;
   case objectValue:
      return value_.map_->compare( *(other.value_.map_) ) < 0;
#endif
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return false;  // unreachable
}

bool
Value::operator <=( const Value &other ) const {
   return !(other < *this);
}

bool
Value::operator >=( const Value &other ) const {
   return !(*this < other);
}

bool
Value::operator >( const Value &other ) const {
   return other < *this;
}

bool
Value::operator ==( const Value &other ) const {
   //if ( type_ != other.type_ )
   // GCC 2.95.3 says:
   // attempt to take address of bit-field structure member `Json::Value::type_'
   // Beats me, but a temp solves the problem.
   int temp = other.type_;
   if ( type_ != temp )
      return false;
   switch ( type_ ) {
   case nullValue:
      return true;
   case intValue:
      return value_.int_ == other.value_.int_;
   case uintValue:
      return value_.uint_ == other.value_.uint_;
   case realValue:
      return value_.real_ == other.value_.real_;
   case booleanValue:
      return value_.bool_ == other.value_.bool_;
   case stringValue:
      return ( value_.string_ == other.value_.string_ )
             || ( other.value_.string_
                  &&  value_.string_
                  && strcmp( value_.string_, other.value_.string_ ) == 0 );
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
   case objectValue:
      return value_.map_->size() == other.value_.map_->size()
             && (*value_.map_) == (*other.value_.map_);
#else
   case arrayValue:
      return value_.array_->compare( *(other.value_.array_) ) == 0;
   case objectValue:
      return value_.map_->compare( *(other.value_.map_) ) == 0;
#endif
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return false;  // unreachable
}

bool
Value::operator !=( const Value &other ) const
{
   return !( *this == other );
}

const char *
Value::asCString() const {
   JSON_ASSERT( type_ == stringValue );
   return value_.string_;
}


std::string
Value::asString() const {
   switch ( type_ ) {
   case nullValue:
      return "";
   case stringValue:
      return value_.string_ ? value_.string_ : "";
   case booleanValue:
      return value_.bool_ ? "true" : "false";
   case intValue:
   case uintValue:
   case realValue:
   case arrayValue:
   case objectValue:
      JSON_FAIL_MESSAGE( "Type is not convertible to string" );
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return ""; // unreachable
}

# ifdef JSON_USE_CPPTL
CppTL::ConstString
Value::asConstString() const
{
   return CppTL::ConstString( asString().c_str() );
}
# endif


Value::Int
Value::asInt() const {
   switch ( type_ )
   {
   case nullValue:
      return 0;
   case intValue:
      JSON_ASSERT_MESSAGE( value_.int_ >= minInt  &&  value_.int_ <= maxInt, "unsigned integer out of signed int range" );
      return Int(value_.int_);
   case uintValue:
      JSON_ASSERT_MESSAGE( value_.uint_ <= UInt(maxInt), "unsigned integer out of signed int range" );
      return Int(value_.uint_);
   case realValue:
      JSON_ASSERT_MESSAGE( value_.real_ >= minInt  &&  value_.real_ <= maxInt, "Real out of signed integer range" );
      return Int( value_.real_ );
   case booleanValue:
      return value_.bool_ ? 1 : 0;
   case stringValue:
   case arrayValue:
   case objectValue:
      JSON_FAIL_MESSAGE( "Type is not convertible to int" );
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return 0; // unreachable;
}


Value::UInt
Value::asUInt() const {
   switch ( type_ ) {
   case nullValue:
      return 0;
   case intValue:
      JSON_ASSERT_MESSAGE( value_.int_ >= 0, "Negative integer can not be converted to unsigned integer" );
      JSON_ASSERT_MESSAGE( value_.int_ <= maxUInt, "signed integer out of UInt range" );
      return UInt(value_.int_);
   case uintValue:
      JSON_ASSERT_MESSAGE( value_.uint_ <= maxUInt, "unsigned integer out of UInt range" );
      return UInt(value_.uint_);
   case realValue:
      JSON_ASSERT_MESSAGE( value_.real_ >= 0  &&  value_.real_ <= maxUInt,  "Real out of unsigned integer range" );
      return UInt( value_.real_ );
   case booleanValue:
      return value_.bool_ ? 1 : 0;
   case stringValue:
   case arrayValue:
   case objectValue:
      JSON_FAIL_MESSAGE( "Type is not convertible to uint" );
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return 0; // unreachable;
}


# if defined(JSON_HAS_INT64)

Value::Int64
Value::asInt64() const {
   switch ( type_ ) {
   case nullValue:
      return 0;
   case intValue:
      return value_.int_;
   case uintValue:
      JSON_ASSERT_MESSAGE( value_.uint_ <= UInt64(maxInt64), "unsigned integer out of Int64 range" );
      return value_.uint_;
   case realValue:
      JSON_ASSERT_MESSAGE( value_.real_ >= minInt64  &&  value_.real_ <= maxInt64, "Real out of Int64 range" );
      return Int( value_.real_ );
   case booleanValue:
      return value_.bool_ ? 1 : 0;
   case stringValue:
   case arrayValue:
   case objectValue:
      JSON_FAIL_MESSAGE( "Type is not convertible to Int64" );
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return 0; // unreachable;
}


Value::UInt64
Value::asUInt64() const {
   switch ( type_ ) {
   case nullValue:
      return 0;
   case intValue:
      JSON_ASSERT_MESSAGE( value_.int_ >= 0, "Negative integer can not be converted to UInt64" );
      return value_.int_;
   case uintValue:
      return value_.uint_;
   case realValue:
      JSON_ASSERT_MESSAGE( value_.real_ >= 0  &&  value_.real_ <= maxUInt64,  "Real out of UInt64 range" );
      return UInt( value_.real_ );
   case booleanValue:
      return value_.bool_ ? 1 : 0;
   case stringValue:
   case arrayValue:
   case objectValue:
      JSON_FAIL_MESSAGE( "Type is not convertible to UInt64" );
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return 0; // unreachable;
}
# endif // if defined(JSON_HAS_INT64)


LargestInt
Value::asLargestInt() const {
#if defined(JSON_NO_INT64)
    return asInt();
#else
    return asInt64();
#endif
}


LargestUInt
Value::asLargestUInt() const {
#if defined(JSON_NO_INT64)
    return asUInt();
#else
    return asUInt64();
#endif
}


double
Value::asDouble() const {
   switch ( type_ ) {
   case nullValue:
      return 0.0;
   case intValue:
      return static_cast<double>( value_.int_ );
   case uintValue:
#if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
      return static_cast<double>( value_.uint_ );
#else // if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
      return static_cast<double>( Int(value_.uint_/2) ) * 2 + Int(value_.uint_ & 1);
#endif // if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
   case realValue:
      return value_.real_;
   case booleanValue:
      return value_.bool_ ? 1.0 : 0.0;
   case stringValue:
   case arrayValue:
   case objectValue:
      JSON_FAIL_MESSAGE( "Type is not convertible to double" );
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return 0; // unreachable;
}

float
Value::asFloat() const {
   switch ( type_ ) {
   case nullValue:
      return 0.0f;
   case intValue:
      return static_cast<float>( value_.int_ );
   case uintValue:
#if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
      return static_cast<float>( value_.uint_ );
#else // if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
      return static_cast<float>( Int(value_.uint_/2) ) * 2 + Int(value_.uint_ & 1);
#endif // if !defined(JSON_USE_INT64_DOUBLE_CONVERSION)
   case realValue:
      return static_cast<float>( value_.real_ );
   case booleanValue:
      return value_.bool_ ? 1.0f : 0.0f;
   case stringValue:
   case arrayValue:
   case objectValue:
      JSON_FAIL_MESSAGE( "Type is not convertible to float" );
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return 0.0f; // unreachable;
}

bool
Value::asBool() const {
   switch ( type_ ) {
   case nullValue:
      return false;
   case intValue:
   case uintValue:
      return value_.int_ != 0;
   case realValue:
      return value_.real_ != 0.0;
   case booleanValue:
      return value_.bool_;
   case stringValue:
      return value_.string_  &&  value_.string_[0] != 0;
   case arrayValue:
   case objectValue:
      return value_.map_->size() != 0;
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return false; // unreachable;
}


bool
Value::isConvertibleTo( ValueType other ) const {
   switch ( type_ ) {
   case nullValue:
      return true;
   case intValue:
      return ( other == nullValue  &&  value_.int_ == 0 )
             || other == intValue
             || ( other == uintValue  && value_.int_ >= 0 )
             || other == realValue
             || other == stringValue
             || other == booleanValue;
   case uintValue:
      return ( other == nullValue  &&  value_.uint_ == 0 )
             || ( other == intValue  && value_.uint_ <= (unsigned)maxInt )
             || other == uintValue
             || other == realValue
             || other == stringValue
             || other == booleanValue;
   case realValue:
      return ( other == nullValue  &&  value_.real_ == 0.0 )
             || ( other == intValue  &&  value_.real_ >= minInt  &&  value_.real_ <= maxInt )
             || ( other == uintValue  &&  value_.real_ >= 0  &&  value_.real_ <= maxUInt )
             || other == realValue
             || other == stringValue
             || other == booleanValue;
   case booleanValue:
      return ( other == nullValue  &&  value_.bool_ == false )
             || other == intValue
             || other == uintValue
             || other == realValue
             || other == stringValue
             || other == booleanValue;
   case stringValue:
      return other == stringValue
             || ( other == nullValue  &&  (!value_.string_  ||  value_.string_[0] == 0) );
   case arrayValue:
      return other == arrayValue
             ||  ( other == nullValue  &&  value_.map_->size() == 0 );
   case objectValue:
      return other == objectValue
             ||  ( other == nullValue  &&  value_.map_->size() == 0 );
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return false; // unreachable;
}


/// Number of values in array or object
ArrayIndex
Value::size() const {
   switch ( type_ ) {
   case nullValue:
   case intValue:
   case uintValue:
   case realValue:
   case booleanValue:
   case stringValue:
      return 0;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:  // size of the array is highest index + 1
      if ( !value_.map_->empty() ) {
         ObjectValues::const_iterator itLast = value_.map_->end();
         --itLast;
         return (*itLast).first.index()+1;
      }
      return 0;
   case objectValue:
      return ArrayIndex( value_.map_->size() );
#else
   case arrayValue:
      return Int( value_.array_->size() );
   case objectValue:
      return Int( value_.map_->size() );
#endif
   default:
      JSON_ASSERT_UNREACHABLE;
   }
   return 0; // unreachable;
}


bool
Value::empty() const {
   if ( isNull() || isArray() || isObject() )
      return size() == 0u;
   else
      return false;
}


bool
Value::operator!() const {
   return isNull();
}


void
Value::clear() {
   JSON_ASSERT( type_ == nullValue  ||  type_ == arrayValue  || type_ == objectValue );

   switch ( type_ )
   {
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
   case objectValue:
      value_.map_->clear();
      break;
#else
   case arrayValue:
      value_.array_->clear();
      break;
   case objectValue:
      value_.map_->clear();
      break;
#endif
   default:
      break;
   }
}

void
Value::resize( ArrayIndex newSize ) {
   JSON_ASSERT( type_ == nullValue  ||  type_ == arrayValue );
   if ( type_ == nullValue )
      *this = Value( arrayValue );
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   ArrayIndex oldSize = size();
   if ( newSize == 0 )
      clear();
   else if ( newSize > oldSize )
      (*this)[ newSize - 1 ];
   else {
      for ( ArrayIndex index = newSize; index < oldSize; ++index )
      {
         value_.map_->erase( index );
      }
      assert( size() == newSize );
   }
#else
   value_.array_->resize( newSize );
#endif
}


Value &
Value::operator[]( ArrayIndex index ) {
   JSON_ASSERT( type_ == nullValue  ||  type_ == arrayValue );
   if ( type_ == nullValue )
      *this = Value( arrayValue );
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   CZString key( index );
   ObjectValues::iterator it = value_.map_->lower_bound( key );
   if ( it != value_.map_->end()  &&  (*it).first == key )
      return (*it).second;

   ObjectValues::value_type defaultValue( key, null );
   it = value_.map_->insert( it, defaultValue );
   return (*it).second;
#else
   return value_.array_->resolveReference( index );
#endif
}


Value &
Value::operator[]( int index ) {
   JSON_ASSERT( index >= 0 );
   return (*this)[ ArrayIndex(index) ];
}


const Value &
Value::operator[]( ArrayIndex index ) const {
   JSON_ASSERT( type_ == nullValue  ||  type_ == arrayValue );
   if ( type_ == nullValue )
      return null;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   CZString key( index );
   ObjectValues::const_iterator it = value_.map_->find( key );
   if ( it == value_.map_->end() )
      return null;
   return (*it).second;
#else
   Value *value = value_.array_->find( index );
   return value ? *value : null;
#endif
}


const Value &
Value::operator[]( int index ) const {
   JSON_ASSERT( index >= 0 );
   return (*this)[ ArrayIndex(index) ];
}


Value &
Value::operator[]( const char *key ) {
   return resolveReference( key, false );
}


Value &
Value::resolveReference( const char *key,
                         bool isStatic ) {
   JSON_ASSERT( type_ == nullValue  ||  type_ == objectValue );
   if ( type_ == nullValue )
      *this = Value( objectValue );
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   CZString actualKey( key, isStatic ? CZString::noDuplication
                                     : CZString::duplicateOnCopy );
   ObjectValues::iterator it = value_.map_->lower_bound( actualKey );
   if ( it != value_.map_->end()  &&  (*it).first == actualKey )
      return (*it).second;

   ObjectValues::value_type defaultValue( actualKey, null );
   it = value_.map_->insert( it, defaultValue );
   Value &value = (*it).second;
   return value;
#else
   return value_.map_->resolveReference( key, isStatic );
#endif
}


Value
Value::get( ArrayIndex index,
            const Value &defaultValue ) const {
   const Value *value = &((*this)[index]);
   return value == &null ? defaultValue : *value;
}


bool
Value::isValidIndex( ArrayIndex index ) const {
   return index < size();
}



const Value &
Value::operator[]( const char *key ) const {
   JSON_ASSERT( type_ == nullValue  ||  type_ == objectValue );
   if ( type_ == nullValue )
      return null;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   CZString actualKey( key, CZString::noDuplication );
   ObjectValues::const_iterator it = value_.map_->find( actualKey );
   if ( it == value_.map_->end() )
      return null;
   return (*it).second;
#else
   const Value *value = value_.map_->find( key );
   return value ? *value : null;
#endif
}


Value &
Value::operator[]( const std::string &key ) {
   return (*this)[ key.c_str() ];
}


const Value &
Value::operator[]( const std::string &key ) const {
   return (*this)[ key.c_str() ];
}

Value &
Value::operator[]( const StaticString &key ) {
   return resolveReference( key, true );
}


# ifdef JSON_USE_CPPTL
Value &
Value::operator[]( const CppTL::ConstString &key ) {
   return (*this)[ key.c_str() ];
}


const Value &
Value::operator[]( const CppTL::ConstString &key ) const {
   return (*this)[ key.c_str() ];
}
# endif


Value &
Value::append( const Value &value ) {
   return (*this)[size()] = value;
}


Value
Value::get( const char *key,
            const Value &defaultValue ) const {
   const Value *value = &((*this)[key]);
   return value == &null ? defaultValue : *value;
}


Value
Value::get( const std::string &key,
            const Value &defaultValue ) const {
   return get( key.c_str(), defaultValue );
}

Value
Value::removeMember( const char* key ) {
   JSON_ASSERT( type_ == nullValue  ||  type_ == objectValue );
   if ( type_ == nullValue )
      return null;
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   CZString actualKey( key, CZString::noDuplication );
   ObjectValues::iterator it = value_.map_->find( actualKey );
   if ( it == value_.map_->end() )
      return null;
   Value old(it->second);
   value_.map_->erase(it);
   return old;
#else
   Value *value = value_.map_->find( key );
   if (value){
      Value old(*value);
      value_.map_.remove( key );
      return old;
   } else {
      return null;
   }
#endif
}

Value
Value::removeMember( const std::string &key ) {
   return removeMember( key.c_str() );
}

# ifdef JSON_USE_CPPTL
Value
Value::get( const CppTL::ConstString &key,
            const Value &defaultValue ) const {
   return get( key.c_str(), defaultValue );
}
# endif

bool
Value::isMember( const char *key ) const {
   const Value *value = &((*this)[key]);
   return value != &null;
}


bool
Value::isMember( const std::string &key ) const {
   return isMember( key.c_str() );
}


# ifdef JSON_USE_CPPTL
bool
Value::isMember( const CppTL::ConstString &key ) const {
   return isMember( key.c_str() );
}
#endif

Value::Members
Value::getMemberNames() const {
   JSON_ASSERT( type_ == nullValue  ||  type_ == objectValue );
   if ( type_ == nullValue )
       return Value::Members();
   Members members;
   members.reserve( value_.map_->size() );
#ifndef JSON_VALUE_USE_INTERNAL_MAP
   ObjectValues::const_iterator it = value_.map_->begin();
   ObjectValues::const_iterator itEnd = value_.map_->end();
   for ( ; it != itEnd; ++it )
      members.push_back( std::string( (*it).first.c_str() ) );
#else
   ValueInternalMap::IteratorState it;
   ValueInternalMap::IteratorState itEnd;
   value_.map_->makeBeginIterator( it );
   value_.map_->makeEndIterator( itEnd );
   for ( ; !ValueInternalMap::equals( it, itEnd ); ValueInternalMap::increment(it) )
      members.push_back( std::string( ValueInternalMap::key( it ) ) );
#endif
   return members;
}
//
//# ifdef JSON_USE_CPPTL
//EnumMemberNames
//Value::enumMemberNames() const
//{
//   if ( type_ == objectValue )
//   {
//      return CppTL::Enum::any(  CppTL::Enum::transform(
//         CppTL::Enum::keys( *(value_.map_), CppTL::Type<const CZString &>() ),
//         MemberNamesTransform() ) );
//   }
//   return EnumMemberNames();
//}
//
//
//EnumValues
//Value::enumValues() const
//{
//   if ( type_ == objectValue  ||  type_ == arrayValue )
//      return CppTL::Enum::anyValues( *(value_.map_),
//                                     CppTL::Type<const Value &>() );
//   return EnumValues();
//}
//
//# endif


bool
Value::isNull() const {
   return type_ == nullValue;
}


bool
Value::isBool() const {
   return type_ == booleanValue;
}


bool
Value::isInt() const {
   return type_ == intValue;
}


bool
Value::isUInt() const {
   return type_ == uintValue;
}


bool
Value::isIntegral() const {
   return type_ == intValue
          ||  type_ == uintValue
          ||  type_ == booleanValue;
}


bool
Value::isDouble() const {
   return type_ == realValue;
}


bool
Value::isNumeric() const {
   return isIntegral() || isDouble();
}


bool
Value::isString() const {
   return type_ == stringValue;
}


bool
Value::isArray() const {
   return type_ == nullValue  ||  type_ == arrayValue;
}


bool
Value::isObject() const {
   return type_ == nullValue  ||  type_ == objectValue;
}


void
Value::setComment( const char *comment,
                   CommentPlacement placement ) {
   if ( !comments_ )
      comments_ = new CommentInfo[numberOfCommentPlacement];
   comments_[placement].setComment( comment );
}


void
Value::setComment( const std::string &comment,
                   CommentPlacement placement ) {
   setComment( comment.c_str(), placement );
}


bool
Value::hasComment( CommentPlacement placement ) const {
   return comments_ != 0  &&  comments_[placement].comment_ != 0;
}

std::string
Value::getComment( CommentPlacement placement ) const {
   if ( hasComment(placement) )
      return comments_[placement].comment_;
   return "";
}


std::string
Value::toStyledString() const {
   StyledWriter writer;
   return writer.write( *this );
}


Value::const_iterator
Value::begin() const {
   switch ( type_ ) {
#ifdef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
      if ( value_.array_ ) {
         ValueInternalArray::IteratorState it;
         value_.array_->makeBeginIterator( it );
         return const_iterator( it );
      }
      break;
   case objectValue:
      if ( value_.map_ ) {
         ValueInternalMap::IteratorState it;
         value_.map_->makeBeginIterator( it );
         return const_iterator( it );
      }
      break;
#else
   case arrayValue:
   case objectValue:
      if ( value_.map_ )
         return const_iterator( value_.map_->begin() );
      break;
#endif
   default:
      break;
   }
   return const_iterator();
}

Value::const_iterator
Value::end() const {
   switch ( type_ ) {
#ifdef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
      if ( value_.array_ ) {
         ValueInternalArray::IteratorState it;
         value_.array_->makeEndIterator( it );
         return const_iterator( it );
      }
      break;
   case objectValue:
      if ( value_.map_ ) {
         ValueInternalMap::IteratorState it;
         value_.map_->makeEndIterator( it );
         return const_iterator( it );
      }
      break;
#else
   case arrayValue:
   case objectValue:
      if ( value_.map_ )
         return const_iterator( value_.map_->end() );
      break;
#endif
   default:
      break;
   }
   return const_iterator();
}


Value::iterator
Value::begin() {
   switch ( type_ ) {
#ifdef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
      if ( value_.array_ ) {
         ValueInternalArray::IteratorState it;
         value_.array_->makeBeginIterator( it );
         return iterator( it );
      }
      break;
   case objectValue:
      if ( value_.map_ ) {
         ValueInternalMap::IteratorState it;
         value_.map_->makeBeginIterator( it );
         return iterator( it );
      }
      break;
#else
   case arrayValue:
   case objectValue:
      if ( value_.map_ )
         return iterator( value_.map_->begin() );
      break;
#endif
   default:
      break;
   }
   return iterator();
}

Value::iterator
Value::end() {
   switch ( type_ ) {
#ifdef JSON_VALUE_USE_INTERNAL_MAP
   case arrayValue:
      if ( value_.array_ ) {
         ValueInternalArray::IteratorState it;
         value_.array_->makeEndIterator( it );
         return iterator( it );
      }
      break;
   case objectValue:
      if ( value_.map_ ) {
         ValueInternalMap::IteratorState it;
         value_.map_->makeEndIterator( it );
         return iterator( it );
      }
      break;
#else
   case arrayValue:
   case objectValue:
      if ( value_.map_ )
         return iterator( value_.map_->end() );
      break;
#endif
   default:
      break;
   }
   return iterator();
}


// class PathArgument
// //////////////////////////////////////////////////////////////////

PathArgument::PathArgument()
   : kind_( kindNone ) {
}


PathArgument::PathArgument( ArrayIndex index )
   : index_( index )
   , kind_( kindIndex ) {
}


PathArgument::PathArgument( const char *key )
   : key_( key )
   , kind_( kindKey ) {
}


PathArgument::PathArgument( const std::string &key )
   : key_( key.c_str() )
   , kind_( kindKey ) {
}

// class Path
// //////////////////////////////////////////////////////////////////

Path::Path( const std::string &path,
            const PathArgument &a1,
            const PathArgument &a2,
            const PathArgument &a3,
            const PathArgument &a4,
            const PathArgument &a5 ) {
   InArgs in;
   in.push_back( &a1 );
   in.push_back( &a2 );
   in.push_back( &a3 );
   in.push_back( &a4 );
   in.push_back( &a5 );
   makePath( path, in );
}


void
Path::makePath( const std::string &path,
                const InArgs &in ) {
   const char *current = path.c_str();
   const char *end = current + path.length();
   InArgs::const_iterator itInArg = in.begin();
   while ( current != end ) {
      if ( *current == '[' ) {
         ++current;
         if ( *current == '%' )
            addPathInArg( path, in, itInArg, PathArgument::kindIndex );
         else {
            ArrayIndex index = 0;
            for ( ; current != end && *current >= '0'  &&  *current <= '9'; ++current )
               index = index * 10 + ArrayIndex(*current - '0');
            args_.push_back( index );
         }
         if ( current == end  ||  *current++ != ']' )
            invalidPath( path, int(current - path.c_str()) );
      }
      else if ( *current == '%' ) {
         addPathInArg( path, in, itInArg, PathArgument::kindKey );
         ++current;
      }
      else if ( *current == '.' ) {
         ++current;
      } else {
         const char *beginName = current;
         while ( current != end  &&  !strchr( "[.", *current ) )
            ++current;
         args_.push_back( std::string( beginName, current ) );
      }
   }
}


void
Path::addPathInArg( const std::string &path,
                    const InArgs &in,
                    InArgs::const_iterator &itInArg,
                    PathArgument::Kind kind ) {
   if ( itInArg == in.end() ) {
      // Error: missing argument %d
   } else if ( (*itInArg)->kind_ != kind ) {
      // Error: bad argument type
   } else {
      args_.push_back( **itInArg );
   }
}


void
Path::invalidPath( const std::string &path,
                   int location ) {
   // Error: invalid path.
}


const Value &
Path::resolve( const Value &root ) const {
   const Value *node = &root;
   for ( Args::const_iterator it = args_.begin(); it != args_.end(); ++it ) {
      const PathArgument &arg = *it;
      if ( arg.kind_ == PathArgument::kindIndex ) {
         if ( !node->isArray()  ||  node->isValidIndex( arg.index_ ) ) {
            // Error: unable to resolve path (array value expected at position...
         }
         node = &((*node)[arg.index_]);
      } else if ( arg.kind_ == PathArgument::kindKey ) {
         if ( !node->isObject() ) {
            // Error: unable to resolve path (object value expected at position...)
         }
         node = &((*node)[arg.key_]);
         if ( node == &Value::null ) {
            // Error: unable to resolve path (object has no member named '' at position...)
         }
      }
   }
   return *node;
}


Value
Path::resolve( const Value &root,
               const Value &defaultValue ) const {
   const Value *node = &root;
   for ( Args::const_iterator it = args_.begin(); it != args_.end(); ++it ) {
      const PathArgument &arg = *it;
      if ( arg.kind_ == PathArgument::kindIndex ) {
         if ( !node->isArray()  ||  node->isValidIndex( arg.index_ ) )
            return defaultValue;
         node = &((*node)[arg.index_]);
      } else if ( arg.kind_ == PathArgument::kindKey ) {
         if ( !node->isObject() )
            return defaultValue;
         node = &((*node)[arg.key_]);
         if ( node == &Value::null )
            return defaultValue;
      }
   }
   return *node;
}


Value &
Path::make( Value &root ) const {
   Value *node = &root;
   for ( Args::const_iterator it = args_.begin(); it != args_.end(); ++it ) {
      const PathArgument &arg = *it;
      if ( arg.kind_ == PathArgument::kindIndex ) {
         if ( !node->isArray() ) {
            // Error: node is not an array at position ...
         }
         node = &((*node)[arg.index_]);
      } else if ( arg.kind_ == PathArgument::kindKey ) {
         if ( !node->isObject() ) {
            // Error: node is not an object at position...
         }
         node = &((*node)[arg.key_]);
      }
   }
   return *node;
}


} // namespace Json

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_value.cpp
// //////////////////////////////////////////////////////////////////////






// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: src/lib_json/json_writer.cpp
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#if !defined(JSON_IS_AMALGAMATION)
# include <json/writer.h>
# include "json_tool.h"
#endif // if !defined(JSON_IS_AMALGAMATION)
#include <utility>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#if _MSC_VER >= 1400 // VC++ 8.0
#pragma warning( disable : 4996 )   // disable warning about strdup being deprecated.
#endif

namespace Json {

static bool containsControlCharacter( const char* str ) {
   while ( *str ) {
      if ( isControlCharacter( *(str++) ) )
         return true;
   }
   return false;
}


std::string valueToString( LargestInt value ) {
   UIntToStringBuffer buffer;
   char *current = buffer + sizeof(buffer);
   bool isNegative = value < 0;
   if ( isNegative )
      value = -value;
   uintToString( LargestUInt(value), current );
   if ( isNegative )
      *--current = '-';
   assert( current >= buffer );
   return current;
}


std::string valueToString( LargestUInt value ) {
   UIntToStringBuffer buffer;
   char *current = buffer + sizeof(buffer);
   uintToString( value, current );
   assert( current >= buffer );
   return current;
}

#if defined(JSON_HAS_INT64)

std::string valueToString( Int value )
{
   return valueToString( LargestInt(value) );
}


std::string valueToString( UInt value ) {
   return valueToString( LargestUInt(value) );
}

#endif // # if defined(JSON_HAS_INT64)


std::string valueToString( double value ) {
   char buffer[32];
#if defined(_MSC_VER) && defined(__STDC_SECURE_LIB__) // Use secure version with visual studio 2005 to avoid warning.
   sprintf_s(buffer, sizeof(buffer), "%#.16g", value);
#else
   sprintf(buffer, "%#.16g", value);
#endif
   char* ch = buffer + strlen(buffer) - 1;
   if (*ch != '0') return buffer; // nothing to truncate, so save time
   while(ch > buffer && *ch == '0'){
     --ch;
   }
   char* last_nonzero = ch;
   while(ch >= buffer){
     switch(*ch){
     case '0':
     case '1':
     case '2':
     case '3':
     case '4':
     case '5':
     case '6':
     case '7':
     case '8':
     case '9':
       --ch;
       continue;
     case '.':
       // Truncate zeroes to save bytes in output, but keep one.
       *(last_nonzero+2) = '\0';
       return buffer;
     default:
       return buffer;
     }
   }
   return buffer;
}


std::string valueToString( bool value ) {
   return value ? "true" : "false";
}

std::string valueToQuotedString( const char *value ) {
   // Not sure how to handle unicode...
   if (strpbrk(value, "\"\\\b\f\n\r\t") == NULL && !containsControlCharacter( value ))
      return std::string("\"") + value + "\"";
   // We have to walk value and escape any special characters.
   // Appending to std::string is not efficient, but this should be rare.
   // (Note: forward slashes are *not* rare, but I am not escaping them.)
   std::string::size_type maxsize = strlen(value)*2 + 3; // allescaped+quotes+NULL
   std::string result;
   result.reserve(maxsize); // to avoid lots of mallocs
   result += "\"";
   for (const char* c=value; *c != 0; ++c) {
      switch(*c) {
         case '\"':
            result += "\\\"";
            break;
         case '\\':
            result += "\\\\";
            break;
         case '\b':
            result += "\\b";
            break;
         case '\f':
            result += "\\f";
            break;
         case '\n':
            result += "\\n";
            break;
         case '\r':
            result += "\\r";
            break;
         case '\t':
            result += "\\t";
            break;
         //case '/':
            // Even though \/ is considered a legal escape in JSON, a bare
            // slash is also legal, so I see no reason to escape it.
            // (I hope I am not misunderstanding something.
            // blep notes: actually escaping \/ may be useful in javascript to avoid </
            // sequence.
            // Should add a flag to allow this compatibility mode and prevent this
            // sequence from occurring.
         default:
            if ( isControlCharacter( *c ) ) {
               std::ostringstream oss;
               oss << "\\u" << std::hex << std::uppercase << std::setfill('0') << std::setw(4) << static_cast<int>(*c);
               result += oss.str();
            } else {
               result += *c;
            }
            break;
      }
   }
   result += "\"";
   return result;
}

// Class Writer
// //////////////////////////////////////////////////////////////////
Writer::~Writer() {
}


// Class FastWriter
// //////////////////////////////////////////////////////////////////

FastWriter::FastWriter()
   : yamlCompatiblityEnabled_( false ) {
}


void
FastWriter::enableYAMLCompatibility() {
   yamlCompatiblityEnabled_ = true;
}


std::string
FastWriter::write( const Value &root ) {
   document_ = "";
   writeValue( root );
   document_ += "\n";
   return document_;
}


void
FastWriter::writeValue( const Value &value ) {
   switch ( value.type() ) {
   case nullValue:
      document_ += "null";
      break;
   case intValue:
      document_ += valueToString( value.asLargestInt() );
      break;
   case uintValue:
      document_ += valueToString( value.asLargestUInt() );
      break;
   case realValue:
      document_ += valueToString( value.asDouble() );
      break;
   case stringValue:
      document_ += valueToQuotedString( value.asCString() );
      break;
   case booleanValue:
      document_ += valueToString( value.asBool() );
      break;
   case arrayValue: {
         document_ += "[";
         int size = value.size();
         for ( int index =0; index < size; ++index ) {
            if ( index > 0 )
               document_ += ",";
            writeValue( value[index] );
         }
         document_ += "]";
      }
      break;
   case objectValue: {
         Value::Members members( value.getMemberNames() );
         document_ += "{";
         for ( Value::Members::iterator it = members.begin();
               it != members.end();
               ++it ) {
            const std::string &name = *it;
            if ( it != members.begin() )
               document_ += ",";
            document_ += valueToQuotedString( name.c_str() );
            document_ += yamlCompatiblityEnabled_ ? ": "
                                                  : ":";
            writeValue( value[name] );
         }
         document_ += "}";
      }
      break;
   }
}


// Class StyledWriter
// //////////////////////////////////////////////////////////////////

StyledWriter::StyledWriter()
   : rightMargin_( 74 )
   , indentSize_( 3 ) {
}


std::string
StyledWriter::write( const Value &root ) {
   document_ = "";
   addChildValues_ = false;
   indentString_ = "";
   writeCommentBeforeValue( root );
   writeValue( root );
   writeCommentAfterValueOnSameLine( root );
   document_ += "\n";
   return document_;
}


void
StyledWriter::writeValue( const Value &value ) {
   switch ( value.type() ) {
   case nullValue:
      pushValue( "null" );
      break;
   case intValue:
      pushValue( valueToString( value.asLargestInt() ) );
      break;
   case uintValue:
      pushValue( valueToString( value.asLargestUInt() ) );
      break;
   case realValue:
      pushValue( valueToString( value.asDouble() ) );
      break;
   case stringValue:
      pushValue( valueToQuotedString( value.asCString() ) );
      break;
   case booleanValue:
      pushValue( valueToString( value.asBool() ) );
      break;
   case arrayValue:
      writeArrayValue( value);
      break;
   case objectValue:
      {
         Value::Members members( value.getMemberNames() );
         if ( members.empty() )
            pushValue( "{}" );
         else {
            writeWithIndent( "{" );
            indent();
            Value::Members::iterator it = members.begin();
            for (;;) {
               const std::string &name = *it;
               const Value &childValue = value[name];
               writeCommentBeforeValue( childValue );
               writeWithIndent( valueToQuotedString( name.c_str() ) );
               document_ += " : ";
               writeValue( childValue );
               if ( ++it == members.end() ) {
                  writeCommentAfterValueOnSameLine( childValue );
                  break;
               }
               document_ += ",";
               writeCommentAfterValueOnSameLine( childValue );
            }
            unindent();
            writeWithIndent( "}" );
         }
      }
      break;
   }
}


void
StyledWriter::writeArrayValue( const Value &value ) {
   unsigned size = value.size();
   if ( size == 0 )
      pushValue( "[]" );
   else {
      bool isArrayMultiLine = isMultineArray( value );
      if ( isArrayMultiLine ) {
         writeWithIndent( "[" );
         indent();
         bool hasChildValue = !childValues_.empty();
         unsigned index =0;
         for (;;) {
            const Value &childValue = value[index];
            writeCommentBeforeValue( childValue );
            if ( hasChildValue )
               writeWithIndent( childValues_[index] );
            else {
               writeIndent();
               writeValue( childValue );
            }
            if ( ++index == size ) {
               writeCommentAfterValueOnSameLine( childValue );
               break;
            }
            document_ += ",";
            writeCommentAfterValueOnSameLine( childValue );
         }
         unindent();
         writeWithIndent( "]" );
      } else { // output on a single line
         assert( childValues_.size() == size );
         document_ += "[ ";
         for ( unsigned index =0; index < size; ++index ) {
            if ( index > 0 )
               document_ += ", ";
            document_ += childValues_[index];
         }
         document_ += " ]";
      }
   }
}


bool
StyledWriter::isMultineArray( const Value &value ) {
   int size = value.size();
   bool isMultiLine = size*3 >= rightMargin_ ;
   childValues_.clear();
   for ( int index =0; index < size  &&  !isMultiLine; ++index ) {
      const Value &childValue = value[index];
      isMultiLine = isMultiLine  ||
                     ( (childValue.isArray()  ||  childValue.isObject())  &&
                        childValue.size() > 0 );
   }
   if ( !isMultiLine ) { // check if line length > max line length
      childValues_.reserve( size );
      addChildValues_ = true;
      int lineLength = 4 + (size-1)*2; // '[ ' + ', '*n + ' ]'
      for ( int index =0; index < size  &&  !isMultiLine; ++index ) {
         writeValue( value[index] );
         lineLength += int( childValues_[index].length() );
         isMultiLine = isMultiLine  &&  hasCommentForValue( value[index] );
      }
      addChildValues_ = false;
      isMultiLine = isMultiLine  ||  lineLength >= rightMargin_;
   }
   return isMultiLine;
}


void
StyledWriter::pushValue( const std::string &value ) {
   if ( addChildValues_ )
      childValues_.push_back( value );
   else
      document_ += value;
}


void
StyledWriter::writeIndent() {
   if ( !document_.empty() )
   {
      char last = document_[document_.length()-1];
      if ( last == ' ' )     // already indented
         return;
      if ( last != '\n' )    // Comments may add new-line
         document_ += '\n';
   }
   document_ += indentString_;
}


void
StyledWriter::writeWithIndent( const std::string &value ) {
   writeIndent();
   document_ += value;
}


void
StyledWriter::indent() {
   indentString_ += std::string( indentSize_, ' ' );
}


void
StyledWriter::unindent() {
   assert( int(indentString_.size()) >= indentSize_ );
   indentString_.resize( indentString_.size() - indentSize_ );
}


void
StyledWriter::writeCommentBeforeValue( const Value &root ) {
   if ( !root.hasComment( commentBefore ) )
      return;
   document_ += normalizeEOL( root.getComment( commentBefore ) );
   document_ += "\n";
}


void
StyledWriter::writeCommentAfterValueOnSameLine( const Value &root ) {
   if ( root.hasComment( commentAfterOnSameLine ) )
      document_ += " " + normalizeEOL( root.getComment( commentAfterOnSameLine ) );

   if ( root.hasComment( commentAfter ) ) {
      document_ += "\n";
      document_ += normalizeEOL( root.getComment( commentAfter ) );
      document_ += "\n";
   }
}


bool
StyledWriter::hasCommentForValue( const Value &value ) {
   return value.hasComment( commentBefore )
          ||  value.hasComment( commentAfterOnSameLine )
          ||  value.hasComment( commentAfter );
}


std::string
StyledWriter::normalizeEOL( const std::string &text ) {
   std::string normalized;
   normalized.reserve( text.length() );
   const char *begin = text.c_str();
   const char *end = begin + text.length();
   const char *current = begin;
   while ( current != end ) {
      char c = *current++;
      if ( c == '\r' ) { // mac or dos EOL
         if ( *current == '\n' ) // convert dos EOL
            ++current;
         normalized += '\n';
      } else // handle unix EOL & other char
         normalized += c;
   }
   return normalized;
}


// Class StyledStreamWriter
// //////////////////////////////////////////////////////////////////

StyledStreamWriter::StyledStreamWriter( std::string indentation )
   : document_(NULL)
   , rightMargin_( 74 )
   , indentation_( indentation ) {
}


void
StyledStreamWriter::write( std::ostream &out, const Value &root ) {
   document_ = &out;
   addChildValues_ = false;
   indentString_ = "";
   writeCommentBeforeValue( root );
   writeValue( root );
   writeCommentAfterValueOnSameLine( root );
   *document_ << "\n";
   document_ = NULL; // Forget the stream, for safety.
}


void
StyledStreamWriter::writeValue( const Value &value ) {
   switch ( value.type() ) {
   case nullValue:
      pushValue( "null" );
      break;
   case intValue:
      pushValue( valueToString( value.asLargestInt() ) );
      break;
   case uintValue:
      pushValue( valueToString( value.asLargestUInt() ) );
      break;
   case realValue:
      pushValue( valueToString( value.asDouble() ) );
      break;
   case stringValue:
      pushValue( valueToQuotedString( value.asCString() ) );
      break;
   case booleanValue:
      pushValue( valueToString( value.asBool() ) );
      break;
   case arrayValue:
      writeArrayValue( value);
      break;
   case objectValue: {
         Value::Members members( value.getMemberNames() );
         if ( members.empty() )
            pushValue( "{}" );
         else {
            writeWithIndent( "{" );
            indent();
            Value::Members::iterator it = members.begin();
            for (;;) {
               const std::string &name = *it;
               const Value &childValue = value[name];
               writeCommentBeforeValue( childValue );
               writeWithIndent( valueToQuotedString( name.c_str() ) );
               *document_ << " : ";
               writeValue( childValue );
               if ( ++it == members.end() ) {
                  writeCommentAfterValueOnSameLine( childValue );
                  break;
               }
               *document_ << ",";
               writeCommentAfterValueOnSameLine( childValue );
            }
            unindent();
            writeWithIndent( "}" );
         }
      }
      break;
   }
}


void
StyledStreamWriter::writeArrayValue( const Value &value ) {
   unsigned size = value.size();
   if ( size == 0 )
      pushValue( "[]" );
   else {
      bool isArrayMultiLine = isMultineArray( value );
      if ( isArrayMultiLine ) {
         writeWithIndent( "[" );
         indent();
         bool hasChildValue = !childValues_.empty();
         unsigned index =0;
         for (;;) {
            const Value &childValue = value[index];
            writeCommentBeforeValue( childValue );
            if ( hasChildValue )
               writeWithIndent( childValues_[index] );
            else {
               writeIndent();
               writeValue( childValue );
            }
            if ( ++index == size ) {
               writeCommentAfterValueOnSameLine( childValue );
               break;
            }
            *document_ << ",";
            writeCommentAfterValueOnSameLine( childValue );
         }
         unindent();
         writeWithIndent( "]" );
      }
      else { // output on a single line
         assert( childValues_.size() == size );
         *document_ << "[ ";
         for ( unsigned index =0; index < size; ++index ) {
            if ( index > 0 )
               *document_ << ", ";
            *document_ << childValues_[index];
         }
         *document_ << " ]";
      }
   }
}


bool
StyledStreamWriter::isMultineArray( const Value &value ) {
   int size = value.size();
   bool isMultiLine = size*3 >= rightMargin_ ;
   childValues_.clear();
   for ( int index =0; index < size  &&  !isMultiLine; ++index ) {
      const Value &childValue = value[index];
      isMultiLine = isMultiLine  ||
                     ( (childValue.isArray()  ||  childValue.isObject())  &&
                        childValue.size() > 0 );
   }
   if ( !isMultiLine ) { // check if line length > max line length
      childValues_.reserve( size );
      addChildValues_ = true;
      int lineLength = 4 + (size-1)*2; // '[ ' + ', '*n + ' ]'
      for ( int index =0; index < size  &&  !isMultiLine; ++index ) {
         writeValue( value[index] );
         lineLength += int( childValues_[index].length() );
         isMultiLine = isMultiLine  &&  hasCommentForValue( value[index] );
      }
      addChildValues_ = false;
      isMultiLine = isMultiLine  ||  lineLength >= rightMargin_;
   }
   return isMultiLine;
}


void
StyledStreamWriter::pushValue( const std::string &value ) {
   if ( addChildValues_ )
      childValues_.push_back( value );
   else
      *document_ << value;
}


void
StyledStreamWriter::writeIndent() {
  /*
    Some comments in this method would have been nice. ;-)

   if ( !document_.empty() )
   {
      char last = document_[document_.length()-1];
      if ( last == ' ' )     // already indented
         return;
      if ( last != '\n' )    // Comments may add new-line
         *document_ << '\n';
   }
  */
   *document_ << '\n' << indentString_;
}


void
StyledStreamWriter::writeWithIndent( const std::string &value ) {
   writeIndent();
   *document_ << value;
}


void
StyledStreamWriter::indent() {
   indentString_ += indentation_;
}


void
StyledStreamWriter::unindent() {
   assert( indentString_.size() >= indentation_.size() );
   indentString_.resize( indentString_.size() - indentation_.size() );
}


void
StyledStreamWriter::writeCommentBeforeValue( const Value &root ) {
   if ( !root.hasComment( commentBefore ) )
      return;
   *document_ << normalizeEOL( root.getComment( commentBefore ) );
   *document_ << "\n";
}


void
StyledStreamWriter::writeCommentAfterValueOnSameLine( const Value &root ) {
   if ( root.hasComment( commentAfterOnSameLine ) )
      *document_ << " " + normalizeEOL( root.getComment( commentAfterOnSameLine ) );

   if ( root.hasComment( commentAfter ) ) {
      *document_ << "\n";
      *document_ << normalizeEOL( root.getComment( commentAfter ) );
      *document_ << "\n";
   }
}


bool
StyledStreamWriter::hasCommentForValue( const Value &value ) {
   return value.hasComment( commentBefore )
          ||  value.hasComment( commentAfterOnSameLine )
          ||  value.hasComment( commentAfter );
}


std::string
StyledStreamWriter::normalizeEOL( const std::string &text ) {
   std::string normalized;
   normalized.reserve( text.length() );
   const char *begin = text.c_str();
   const char *end = begin + text.length();
   const char *current = begin;
   while ( current != end ) {
      char c = *current++;
      if ( c == '\r' ) { // mac or dos EOL
         if ( *current == '\n' ) // convert dos EOL
            ++current;
         normalized += '\n';
      } else // handle unix EOL & other char
         normalized += c;
   }
   return normalized;
}


std::ostream& operator<<( std::ostream &sout, const Value &root ) {
   Json::StyledStreamWriter writer;
   writer.write(sout, root);
   return sout;
}


} // namespace Json

// //////////////////////////////////////////////////////////////////////
// End of content of file: src/lib_json/json_writer.cpp
// //////////////////////////////////////////////////////////////////////





//
// end of /home/mouginot/work/app/pyne/src/jsoncpp.cpp
//


//
// start of /home/mouginot/work/app/pyne/src/jsoncustomwriter.cpp
//
/**********************************************************************
Copyright (c) 2013 by Matt Swain <m.swain@me.com>

The MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

***********************************************************************/

#ifndef PYNE_IS_AMALGAMATED
  #include "json.h"
  #include "jsoncustomwriter.h"
#endif

namespace Json {

CustomWriter::CustomWriter( std::string opencurly,
                            std::string closecurly,
                            std::string opensquare,
                            std::string closesquare,
                            std::string colon,
                            std::string comma,
                            std::string indent,
                            int maxWidth)
   : opencurly_( opencurly )
   , closecurly_( closecurly )
   , opensquare_( opensquare )
   , closesquare_( closesquare )
   , colon_( colon )
   , comma_( comma )
   , indent_( indent )
   , maxWidth_( maxWidth )
{
}


std::string
CustomWriter::write( const Value &root )
{
   document_ = "";
   indentString_ = "";
   writeValue( root, document_, false );
   document_ += "\n";
   return document_;
}


void
CustomWriter::writeValue( const Value &value, std::string &doc, bool forceSingleLine )
{
   switch ( value.type() )
   {
   case nullValue:
      doc += "null";
      break;
   case intValue:
      doc += valueToString( value.asLargestInt() );
      break;
   case uintValue:
      doc += valueToString( value.asLargestUInt() );
      break;
   case realValue:
      doc += valueToString( value.asDouble() );
      break;
   case stringValue:
      doc += valueToQuotedString( value.asCString() );
      break;
   case booleanValue:
      doc += valueToString( value.asBool() );
      break;
   case arrayValue:
      {
         bool isMulti = false;
         if (!forceSingleLine)
         {
            std::string valLine = "";
            writeValue( value, valLine, true);
            if (valLine.length() > maxWidth_)
            {
               isMulti = true;
            }
            else
            {
               doc += valLine;
               break;
            }
         }
         doc += opensquare_;
         if (isMulti)
            indent();
         for ( int index =0; index < value.size(); ++index )
         {
            if (isMulti)
            {
               doc += "\n";
               doc += indentString_;
            }
            writeValue( value[index], doc, false );
            if ( index < value.size()-1 )
               doc += comma_;
         }
         if (isMulti)
         {
            unindent();
            doc += "\n";
            doc += indentString_;
         }
         doc += closesquare_;
      }
      break;
   case objectValue:
      {
         bool isMulti = false;
         if (!forceSingleLine)
         {
            std::string valLine = "";
            writeValue( value, valLine, true);
            if (valLine.length() > maxWidth_)
            {
               isMulti = true;
            }
            else
            {
               doc += valLine;
               break;
            }
         }
         Value::Members members( value.getMemberNames() );
         doc += opencurly_;
         if (isMulti)
            indent();
         for ( Value::Members::iterator it = members.begin();
               it != members.end();
               ++it )
         {
            if (isMulti)
            {
               doc += "\n";
               doc += indentString_;

            }
            const std::string &name = *it;
            doc += valueToQuotedString( name.c_str() );
            doc += colon_;
            writeValue( value[name], doc, forceSingleLine );
            if ( !(it + 1 == members.end()) )
               doc += comma_;
         }
         if (isMulti)
         {
            unindent();
            doc += "\n";
            doc += indentString_;
         }
         doc += closecurly_;
      }
      break;
   }
}


void
CustomWriter::indent()
{
   indentString_ += indent_;
}


void
CustomWriter::unindent()
{
   int idSize = int(indent_.size());
   int idsSize = int(indentString_.size());
   if (idsSize >= idSize)
      indentString_.resize (idsSize - idSize);
}

}
//
// end of /home/mouginot/work/app/pyne/src/jsoncustomwriter.cpp
//


//
// start of /home/mouginot/work/app/pyne/src/material.cpp
//
// Material.cpp
// The very central Material class
// -- Anthony Scopatz

#include <string>
#include <vector>
#include <iomanip>  // std::setprecision
#include <math.h>   // modf
#include <stdexcept>

#ifndef PYNE_IS_AMALGAMATED
#include "transmuters.h"
#include "material.h"
#endif

// h5wrap template
template double h5wrap::get_array_index(hid_t, int, hid_t);
const int mcnp_line_length = 79;


/***************************/
/*** Protected Functions ***/
/***************************/

double pyne::Material::get_comp_sum() {
  // Sums the weights in the composition dictionary
  double sum = 0.0;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    sum = sum + i->second;
  }
  return sum;
}



void pyne::Material::norm_comp() {
  double sum = get_comp_sum();
  if (sum != 1.0 && sum != 0.0) {
    for (comp_iter i = comp.begin(); i != comp.end(); i++)
      i->second = i->second / sum;
  }

  if (mass < 0.0)
    mass = sum;
}






void pyne::Material::_load_comp_protocol0(hid_t db, std::string datapath, int row) {
  hid_t matgroup = H5Gopen2(db, datapath.c_str(), H5P_DEFAULT);
  hid_t nucset;
  double nucvalue;
  ssize_t nuckeylen;
  std::string nuckey;

  // get the number of members in the material group
  H5G_info_t group_info;
  H5Gget_info(matgroup, &group_info);
  hsize_t matG = group_info.nlinks;

  // Iterate over datasets in the group.
  for (int matg = 0; matg < matG; matg++) {
    nuckeylen = 1 + H5Lget_name_by_idx(matgroup, ".", H5_INDEX_NAME, H5_ITER_INC, matg,
                                        NULL, 0, H5P_DEFAULT);
    char * nkey = new char[nuckeylen];
    nuckeylen = H5Lget_name_by_idx(matgroup, ".", H5_INDEX_NAME, H5_ITER_INC, matg,
                                    nkey, nuckeylen, H5P_DEFAULT);
    nuckey = nkey;
    nucset = H5Dopen2(matgroup, nkey, H5P_DEFAULT);
    nucvalue = h5wrap::get_array_index<double>(nucset, row);

    if (nuckey == "Mass" || nuckey == "MASS" || nuckey == "mass")
      mass = nucvalue;
    else
      comp[pyne::nucname::id(nuckey)] = nucvalue;

    H5Dclose(nucset);
    delete[] nkey;
  }

  // Set meta data
  atoms_per_molecule = -1.0;
}



void pyne::Material::_load_comp_protocol1(hid_t db, std::string datapath, int row) {
  std::string nucpath;
  hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  // Grab the nucpath
  hid_t nuc_attr = H5Aopen(data_set, "nucpath", H5P_DEFAULT);
  H5A_info_t nuc_info;
  H5Aget_info(nuc_attr, &nuc_info);
  hsize_t nuc_attr_len = nuc_info.data_size;
  hid_t str_attr = H5Tcopy(H5T_C_S1);
  H5Tset_size(str_attr, nuc_attr_len);
  char * nucpathbuf = new char [nuc_attr_len];
  H5Aread(nuc_attr, str_attr, nucpathbuf);
  nucpath = std::string(nucpathbuf, nuc_attr_len);
  delete[] nucpathbuf;
  H5Tclose(str_attr);
  _load_comp_protocol1(db, datapath, nucpath, row);
}

void pyne::Material::_load_comp_protocol1(hid_t db, std::string datapath, std::string nucpath, int row) {
  bool datapath_exists = h5wrap::path_exists(db, datapath);
  bool nucpath_exists = h5wrap::path_exists(db, nucpath);
  
  if (!nucpath_exists) {
    nucpath = datapath + "_" + nucpath.substr(1);  
    nucpath_exists = h5wrap::path_exists(db, nucpath);
  }

  hid_t data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  hsize_t data_offset[1] = {static_cast<hsize_t>(row)};
  if (row < 0) {
    // Handle negative row indices
    hid_t data_space = H5Dget_space(data_set);
    hsize_t data_dims[1];
    H5Sget_simple_extent_dims(data_space, data_dims, NULL);
    data_offset[0] += data_dims[0];
  }

  // Grab the nuclides
  std::vector<int> nuclides = h5wrap::h5_array_to_cpp_vector_1d<int>(db, nucpath, H5T_NATIVE_INT);
  int nuc_size = nuclides.size();
  hsize_t nuc_dims[1] = {static_cast<hsize_t>(nuc_size)};

  // Get the data hyperslab
  hid_t data_hyperslab = H5Dget_space(data_set);
  hsize_t data_count[1] = {1};
  H5Sselect_hyperslab(data_hyperslab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);

  // Get memory space for writing
  hid_t mem_space = H5Screate_simple(1, data_count, NULL);

  // Get material type
  size_t material_data_size = sizeof(pyne::material_data) + sizeof(double)*(nuc_size-1);
  hid_t desc = H5Tcreate(H5T_COMPOUND, material_data_size);
  hid_t comp_values_array_type = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  H5Tinsert(desc, "mass", HOFFSET(pyne::material_data, mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "density", HOFFSET(pyne::material_data, density),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "atoms_per_molecule", HOFFSET(pyne::material_data, atoms_per_mol),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "comp", HOFFSET(pyne::material_data, comp), comp_values_array_type);

  // make the data array, have to over-allocate
  material_data * mat_data = new material_data [material_data_size];

  // Finally, get data and put in on this instance
  H5Dread(data_set, desc, mem_space, data_hyperslab, H5P_DEFAULT, mat_data);

  mass = (*mat_data).mass;
  density = (*mat_data).density;
  atoms_per_molecule = (*mat_data).atoms_per_mol;
  for (int i = 0; i < nuc_size; i++)
    if((double) (*mat_data).comp[i] != 0) {
      comp[nuclides[i]] = (double) (*mat_data).comp[i];
    }
  delete[] mat_data;
  //
  // Get metadata from associated dataset, if available
  //
  std::string attrpath = datapath + "_metadata";
  bool attrpath_exists = h5wrap::path_exists(db, attrpath);
  if (!attrpath_exists)
    return;

  hid_t metadatapace, attrtype, metadataet, metadatalab, attrmemspace;
  int attrrank;
  hvl_t attrdata [1];

  attrtype = H5Tvlen_create(H5T_NATIVE_CHAR);

  // Get the metadata from the file
  metadataet = H5Dopen2(db, attrpath.c_str(), H5P_DEFAULT);
  metadatalab = H5Dget_space(metadataet);
  H5Sselect_hyperslab(metadatalab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
  attrmemspace = H5Screate_simple(1, data_count, NULL);
  H5Dread(metadataet, attrtype, attrmemspace, metadatalab, H5P_DEFAULT, attrdata);

  // convert to in-memory JSON
  Json::Reader reader;
  reader.parse((char *) attrdata[0].p, (char *) attrdata[0].p+attrdata[0].len, metadata, false);

  // close attr data objects
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(metadataet);
  H5Sclose(metadatapace);
  H5Tclose(attrtype);

  // Close out the HDF5 file
  H5Fclose(db);
}





void pyne::Material::from_hdf5(char * filename, char * datapath, int row, int protocol) {
  std::string fname (filename);
  std::string dpath (datapath);
  from_hdf5(fname, dpath, row, protocol);
}



void pyne::Material::from_hdf5(std::string filename, std::string datapath, int row, int protocol) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(filename);

  //Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);
  // Open the database
  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  bool datapath_exists = h5wrap::path_exists(db, datapath);
  if (!datapath_exists)
    throw h5wrap::PathNotFound(filename, datapath);

  // Clear current content
  comp.clear();

  // Load via various protocols
  if (protocol == 0)
    _load_comp_protocol0(db, datapath, row);
  else if (protocol == 1)
    _load_comp_protocol1(db, datapath, row);
  else
    throw pyne::MaterialProtocolError();

  // Close the database
  status = H5Fclose(db);

  // Renormalize the composition, just to be safe.
  norm_comp();
}

void pyne::Material::from_hdf5(char * filename, char * datapath, char * nucpath, int row, int protocol) {
  std::string fname (filename);
  std::string dpath (datapath);
  std::string npath (nucpath);
  from_hdf5(fname, dpath, nucpath, row, protocol);
}

void pyne::Material::from_hdf5(std::string filename, std::string datapath, std::string nucpath, int row, int protocol) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Check that the file is there
  std::string dpath (datapath);
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (!ish5)
    throw h5wrap::FileNotHDF5(filename);

  //Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);
  // Open the database
  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  bool datapath_exists = h5wrap::path_exists(db, datapath);
  if (!datapath_exists)
    throw h5wrap::PathNotFound(filename, datapath);
  bool nucpath_exists = h5wrap::path_exists(db, nucpath);
  if (!nucpath_exists)
    throw h5wrap::PathNotFound(filename, nucpath);

  // Clear current content
  comp.clear();

  // Load via various protocols
  if (protocol == 0) {
    _load_comp_protocol0(db, datapath, row);
  } else if (protocol == 1) {
    _load_comp_protocol1(db, datapath, nucpath, row);
  } else {
    throw pyne::MaterialProtocolError();
  }
  // Close the database
  status = H5Fclose(db);

  // Renormalize the composition, just to be safe.
  norm_comp();
}



void pyne::Material::write_hdf5(char * filename, char * datapath, char * nucpath, float row, int chunksize) {
  std::string fname (filename);
  std::string groupname (datapath);
  std::string nuclist (nucpath);
  write_hdf5(fname, groupname, nuclist, row, chunksize);
}


void pyne::Material::write_hdf5(std::string filename, std::string datapath,
                                std::string nucpath, float row, int chunksize) {
  int row_num = (int) row;

  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  //Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);
  // Create new/open datafile.
  hid_t db;
  if (pyne::file_exists(filename)) {
    bool ish5 = H5Fis_hdf5(filename.c_str());
    if (!ish5)
      throw h5wrap::FileNotHDF5(filename);
    db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
  }
  else
    db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

  //
  // Read in nuclist if available, write it out if not
  //
  bool nucpath_exists = h5wrap::path_exists(db, nucpath);
  // if no nucpath exists use the nested version...
  if (!nucpath_exists) {
    nucpath = datapath + "_" + nucpath.substr(1);  
    nucpath_exists = h5wrap::path_exists(db, nucpath);
  }

  std::vector<int> nuclides;
  int nuc_size;
  hsize_t nuc_dims[1];

  if (nucpath_exists) {
    nuclides = h5wrap::h5_array_to_cpp_vector_1d<int>(db, nucpath, H5T_NATIVE_INT);
    nuc_size = nuclides.size();
    nuc_dims[0] = nuc_size;
  } else {
    nuclides = std::vector<int>();
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++)
      nuclides.push_back(i->first);
    nuc_size = nuclides.size();

    // Create the data if it doesn't exist
    int nuc_data [nuc_size];
    for (int n = 0; n != nuc_size; n++)
      nuc_data[n] = nuclides[n];
    nuc_dims[0] = nuc_size;
    hid_t nuc_space = H5Screate_simple(1, nuc_dims, NULL);
    hid_t nuc_set = H5Dcreate2(db, nucpath.c_str(), H5T_NATIVE_INT, nuc_space,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(nuc_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nuc_data);
    H5Fflush(db, H5F_SCOPE_GLOBAL);
  }


  //
  // Write out the data itself to the file
  //
  hid_t data_set, data_space, data_hyperslab;
  int data_rank = 1;
  hsize_t data_dims[1] = {1};
  hsize_t data_max_dims[1] = {H5S_UNLIMITED};
  hsize_t data_offset[1] = {0};

  size_t material_data_size = sizeof(pyne::material_data) + sizeof(double)*(nuc_size-1);
  hid_t desc = H5Tcreate(H5T_COMPOUND, material_data_size);
  hid_t comp_values_array_type = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, nuc_dims);

  // make the data table type
  H5Tinsert(desc, "mass", HOFFSET(pyne::material_data, mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "density", HOFFSET(pyne::material_data, density),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "atoms_per_molecule", HOFFSET(pyne::material_data, atoms_per_mol),
            H5T_NATIVE_DOUBLE);
  H5Tinsert(desc, "comp", HOFFSET(pyne::material_data, comp),
            comp_values_array_type);

  material_data * mat_data  = new material_data[material_data_size];
  (*mat_data).mass = mass;
  (*mat_data).density = density;
  (*mat_data).atoms_per_mol = atoms_per_molecule;
  for (int n = 0; n != nuc_size; n++) {
    if (0 < comp.count(nuclides[n]))
      (*mat_data).comp[n] = comp[nuclides[n]];
    else
      (*mat_data).comp[n] = 0.0;
  }

  // get / make the data set
  bool datapath_exists = h5wrap::path_exists(db, datapath);
  if (datapath_exists) {
    data_set = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);
    data_space = H5Dget_space(data_set);
    data_rank = H5Sget_simple_extent_dims(data_space, data_dims, data_max_dims);

    // Determine the row size.
    if (std::signbit(row))
      row_num = data_dims[0] + row;  // careful, row is negative

    if (data_dims[0] <= row_num) {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      H5Dset_extent(data_set, data_dims);
    }

    data_offset[0] = row_num;
  } else {
    // Get full space
    data_space = H5Screate_simple(1, data_dims, data_max_dims);

    // Make data set properties to enable chunking
    hid_t data_set_params = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[1] = {static_cast<hsize_t>(chunksize)};
    H5Pset_chunk(data_set_params, 1, chunk_dims);
    H5Pset_deflate(data_set_params, 1);

    // Create the data set
    data_set = H5Dcreate2(db, datapath.c_str(), desc, data_space, H5P_DEFAULT,
                            data_set_params, H5P_DEFAULT);
    H5Dset_extent(data_set, data_dims);

    // Add attribute pointing to nuc path
    hid_t nuc_attr_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(nuc_attr_type, nucpath.length());
    hid_t nuc_attr_space = H5Screate(H5S_SCALAR);
    hid_t nuc_attr = H5Acreate2(data_set, "nucpath", nuc_attr_type, nuc_attr_space,
                                H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(nuc_attr, nuc_attr_type, nucpath.c_str());
    H5Fflush(db, H5F_SCOPE_GLOBAL);
  }

  // Get the data hyperslab
  data_hyperslab = H5Dget_space(data_set);
  hsize_t data_count[1] = {1};
  H5Sselect_hyperslab(data_hyperslab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);

  // Get a memory space for writing
  hid_t mem_space = H5Screate_simple(1, data_count, data_max_dims);

  // Write the row...
  H5Dwrite(data_set, desc, mem_space, data_hyperslab, H5P_DEFAULT, mat_data);

  // Close out the Dataset
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(data_set);
  H5Sclose(data_space);
  H5Tclose(desc);

  //
  // Write out the metadata to the file
  //
  std::string attrpath = datapath + "_metadata";
  hid_t metadatapace, attrtype, metadataet, metadatalab, attrmemspace;
  int attrrank;

  attrtype = H5Tvlen_create(H5T_NATIVE_CHAR);

  // get / make the data set
  bool attrpath_exists = h5wrap::path_exists(db, attrpath);
  if (attrpath_exists) {
    metadataet = H5Dopen2(db, attrpath.c_str(), H5P_DEFAULT);
    metadatapace = H5Dget_space(metadataet);
    attrrank = H5Sget_simple_extent_dims(metadatapace, data_dims, data_max_dims);

    if (data_dims[0] <= row_num) {
      // row == -0, extend to data set so that we can append, or
      // row_num is larger than current dimension, resize to accomodate.
      data_dims[0] = row_num + 1;
      H5Dset_extent(metadataet, data_dims);
    }

    data_offset[0] = row_num;
  } else {
    hid_t metadataetparams;
    hsize_t attrchunkdims [1];

    // Make data set properties to enable chunking
    metadataetparams = H5Pcreate(H5P_DATASET_CREATE);
    attrchunkdims[0] = chunksize;
    H5Pset_chunk(metadataetparams, 1, attrchunkdims);
    H5Pset_deflate(metadataetparams, 1);

    hvl_t attrfillvalue [1];
    attrfillvalue[0].len = 3;
    attrfillvalue[0].p = (char *) "{}\n";
    H5Pset_fill_value(metadataetparams, attrtype, &attrfillvalue);

    // make dataset
    metadatapace = H5Screate_simple(1, data_dims, data_max_dims);
    metadataet = H5Dcreate2(db, attrpath.c_str(), attrtype, metadatapace,
                         H5P_DEFAULT, metadataetparams, H5P_DEFAULT);
    H5Dset_extent(metadataet, data_dims);
    H5Pclose(metadataetparams);
  }

  // set the attr string
  hvl_t attrdata [1];
  Json::FastWriter writer;
  std::string metadatatr = writer.write(metadata);
  attrdata[0].p = (char *) metadatatr.c_str();
  attrdata[0].len = metadatatr.length();

  // write the attr
  metadatalab = H5Dget_space(metadataet);
  H5Sselect_hyperslab(metadatalab, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
  attrmemspace = H5Screate_simple(1, data_count, data_max_dims);
  H5Dwrite(metadataet, attrtype, attrmemspace, metadatalab, H5P_DEFAULT, attrdata);

  // close attr data objects
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  H5Dclose(metadataet);
  H5Sclose(metadatapace);
  H5Tclose(attrtype);

  // Close out the HDF5 file
  H5Fclose(db);
  // Remember the milk!
  // ...by which I mean to deallocate
  delete[] mat_data;
}

std::string pyne::Material::openmc(std::string frac_type) {
  std::ostringstream oss;

  std::set<int> carbon_set; carbon_set.insert(nucname::id("C"));
  pyne::Material temp_mat = this->expand_elements(carbon_set);

  // vars for consistency
  std::string new_quote = "\"";
  std::string end_quote = "\" ";
  std::string indent = "  ";

  // open the material element
  oss << "<material id=" ;

  // add the mat number
  if (temp_mat.metadata.isMember("mat_number")) {
    int mat_num = temp_mat.metadata["mat_number"].asInt();
    oss << new_quote << mat_num << end_quote;
  }
  // mat numbers are required for openmc
  else {
    throw pyne::ValueError("No material number found in metadata. This is not valid for use in OpenMC.");
    oss << new_quote << "?" << end_quote;
  }

  // add name if specified
  if (temp_mat.metadata.isMember("name")) {
    oss << "name=" << new_quote << temp_mat.metadata["name"].asString() << end_quote;
  }
  // close the material tag
  oss << ">";
  // new line
  oss << std::endl;

  //indent
  oss << indent;

  // specify density
  oss << "<density ";
    // if density is negtaive, report to user
  if (temp_mat.density < 0.0) {
    throw pyne::ValueError("A density < 0.0 was found. This is not valid for use in OpenMC.");
  }
  std::string density_str = std::to_string(temp_mat.density);
  // remove trailing zeros
  density_str.erase( density_str.find_last_not_of('0') + 1, std::string::npos);
  oss << "value=" <<  std::fixed << new_quote << density_str << end_quote;
  oss << "units=" << new_quote << "g/cc" << end_quote << "/>";
  // new line
  oss << std::endl;

  std::map<int, double> fracs;
  std::string frac_attrib;
  if(frac_type == "atom") {
    fracs = temp_mat.to_atom_frac();
    frac_attrib = "ao=";
  }
  else {
    fracs = temp_mat.comp;
    frac_attrib = "wo=";
  }

  // add nuclides
  for(comp_map::iterator f = fracs.begin(); f != fracs.end(); f++) {
    if (f->second == 0.0) { continue; }
    //indent
    oss << "  ";
    // start a new nuclide element
    oss << "<nuclide name=" << new_quote;
    oss << pyne::nucname::openmc(f->first);
    oss << end_quote;
    oss << frac_attrib;
    oss << std::setprecision(4) << std::scientific << new_quote << f->second << end_quote;
    oss << "/>";
    // new line
    oss << std::endl;
  }

  // other OpenMC material properties
  if(temp_mat.metadata.isMember("sab")) {
    oss << indent;
    oss << "<sab name=";
    oss << new_quote << temp_mat.metadata["sab"].asString() << end_quote;
    oss << "/>";
    oss << std::endl;
  }

  if(temp_mat.metadata.isMember("temperature")) {
    oss << indent;
    oss << "<temperature>";
    oss << new_quote << temp_mat.metadata["temperature"].asString() << end_quote;
    oss << "</temperature>";
    oss << std::endl;
  }

  if(temp_mat.metadata.isMember("macroscopic")) {
    oss << indent;
    oss << "<macroscopic name=";
    oss << new_quote << temp_mat.metadata["macroscropic"].asString() << end_quote;
    oss << "/>";
    oss << std::endl;
  }

  if(temp_mat.metadata.isMember("isotropic")) {
    oss << indent;
    oss << "<isotropic>";
    oss << new_quote << temp_mat.metadata["isotropic"].asString() << end_quote;
    oss << "</isotropic>";
    oss << std::endl;
  }

  // close the material node
  oss << "</material>" << std::endl;

  return oss.str();
}

///---------------------------------------------------------------------------//
std::string pyne::Material::get_uwuw_name() {
  // standard uwuw material name is : "mat:<Name of Material>/rho:<density>"
  if (! metadata.isMember("name")) {
    pyne::warning("The material has no name");
    return "";
  }
  std::ostringstream uwuw_name;
  uwuw_name << "mat:";
  uwuw_name << metadata["name"].asString();
  if (density > 0) {
    uwuw_name << "/rho:" << std::setprecision(5) << density;
  } else {
    pyne::warning("No Density defined for this Material");
  }

  return uwuw_name.str();
}

///---------------------------------------------------------------------------//
std::string pyne::Material::mcnp(std::string frac_type) {
  //////////////////// Begin card creation ///////////////////////
  std::ostringstream oss;

  std::string comment_prefix = "C ";
  
  // 'name'
  if (metadata.isMember("name")) {
    oss << "C name: " << metadata["name"].asString() << std::endl;
  }
  // 'density'
  if (density != -1.0) {
     std::stringstream ds;
     ds << std::setprecision(5) << std::fixed << "C density = " << density << std::endl;
     oss << ds.str();
  }
  // 'source'
  if (metadata.isMember("source")) {
     oss << "C source: " << metadata["source"].asString() << std::endl;
  }
  // Metadata comments
  if (metadata.isMember("comments")) {
    std::string comment_string = "comments: " + metadata["comments"].asString();
    oss << pyne::comment_line_wrapping(comment_string, comment_prefix, mcnp_line_length).str();
  }

  // Metadata mat_num
  oss << "m";
  if (metadata.isMember("mat_number")) {
    int mat_num = metadata["mat_number"].asInt();
    oss << mat_num << std::endl;
  } else {
    oss << "?" << std::endl;
  }

  // Set up atom or mass frac map
  std::map<int, double> fracs = get_density_frac(frac_type);
  std::string frac_sign = "";

  // write the frac map
  oss << mcnp_frac(fracs, frac_type);

  return oss.str();
}


///---------------------------------------------------------------------------//
std::string pyne::Material::phits(std::string frac_type) {
  //////////////////// Begin card creation ///////////////////////
  std::ostringstream oss;

  std::string comment_prefix = "C ";
  
  // 'name'
  if (metadata.isMember("name")) {
    oss << "C name: " << metadata["name"].asString() << std::endl;
  }
  // Metadata comments
  if (metadata.isMember("comments")) {
    std::string comment_string = "comments: " + metadata["comments"].asString();
    oss << pyne::comment_line_wrapping(comment_string, comment_prefix, mcnp_line_length).str();
  }

  // Metadata mat_num
  oss << "M[ ";
  if (metadata.isMember("mat_number")) {
    int mat_num = metadata["mat_number"].asInt();
    oss << mat_num;
  } else {
    oss << "?";
  }
  oss << " ]" << std::endl;

  // check for metadata
  std::string keyworkds[6] = {"GAS", "ESTEP", "NLIB", "PLIB", "ELIB", "HLIB"};
  for (auto keyword : keyworkds){
    if (metadata.isMember(keyword)){
      oss << "     "<< keyword << "=" << metadata[keyword].asInt() << std::endl;
    }
  }
  // COND should be "<" or "=" or ">" if present
  if (metadata.isMember("COND")){
    oss << "     COND" << metadata["COND"].asString() << "0" << std::endl;
  }

  // Set up atom or mass frac map
  std::map<int, double> fracs = get_density_frac(frac_type);
  std::string frac_sign = "";

  // write the frac map
  oss << mcnp_frac(fracs, frac_type);

  return oss.str();
}

std::string pyne::Material::mcnp_frac(std::map<int, double> fracs, std::string frac_type){

  std::string frac_sign = "";
  if ("atom" != frac_type) {
    frac_sign = "-";
  }
  
  // iterate through frac map
  // This is an awkward pre-C++11 way to put an int to a string
  std::ostringstream oss;
  for(pyne::comp_iter i = fracs.begin(); i != fracs.end(); ++i) {
    if (i->second > 0.0) {
      // Clear first
      std::stringstream ss;
      std::string nucmcnp;
      std::string table_item;
      ss << pyne::nucname::mcnp(i->first);
      nucmcnp = ss.str();

      int mcnp_id;
      mcnp_id = pyne::nucname::mcnp(i->first);
      // Spaces are important for tests
      table_item = metadata["table_ids"][nucmcnp].asString();
      if (!table_item.empty()) {
        oss << "     " << mcnp_id << "." << table_item << " ";
      } else {
        oss << "     " << mcnp_id << " ";
      }
      // The int needs a little formatting
      std::stringstream fs;
      fs << std::setprecision(4) << std::scientific << frac_sign << i->second << std::endl;
      oss << fs.str();
    }
  }
  return oss.str();
}

///---------------------------------------------------------------------------//
/// Create a set out of the static string array.
std::set<std::string> fluka_builtin(pyne::fluka_mat_strings,
                                    pyne::fluka_mat_strings+pyne::FLUKA_MAT_NUM);

///---------------------------------------------------------------------------//
/// not_fluka_builtin
///---------------------------------------------------------------------------//
/// Convenience function
/// This is written as a negative because that is what we care about
bool pyne::Material::not_fluka_builtin(std::string fluka_name) {
  return (fluka_builtin.find(fluka_name) == fluka_builtin.end());
}

///---------------------------------------------------------------------------//
/// fluka
///---------------------------------------------------------------------------//
/// Main external call
std::string pyne::Material::fluka(int id, std::string frac_type) {
  std::stringstream rs;

  // Element, one nucid
  if (comp.size() == 1) {
    rs << fluka_material_str(id);
  } else if (comp.size() > 1) {
  // Compound
    rs << fluka_compound_str(id, frac_type);
  } else {
    rs << "There is no nuclide information in the Material Object" << std::endl;
  }
  return rs.str();
}

///---------------------------------------------------------------------------//
/// fluka_material_str
///---------------------------------------------------------------------------//
///
/// Requirement:  the material upon which this function is called has
///               exactly one nucid component, i.e. it is elemental
/// Do not assume fluka_name is defined in the metadata.  This function
/// may be called from a user-defined material, i.e. on that is not
/// read out of a UW^2-tagged geometry file, and thus does not have
/// certain metadata.
std::string pyne::Material::fluka_material_str(int id) {
  std::stringstream ms;
  std::string fluka_name; // needed to determine if built-in

  int nucid = comp.begin()->first;

  // NOTE:  first part of 'if' may never be called
  if (metadata.isMember("fluka_name")) {
    fluka_name = metadata["fluka_name"].asString();
  } else {  // Should be elemental
    if (comp.size() > 1 ) {
      std::cerr << "Error: this mix is a compound, there should be a fluka_name defined."
                << std::endl;
      return ms.str();
    }
    fluka_name = nucname::fluka(nucid);
  }

  if (not_fluka_builtin(fluka_name)) {
    ms << fluka_material_component(id, nucid, fluka_name);
  }

  // could be empty
  return ms.str();
}

///---------------------------------------------------------------------------//
/// fluka_material_component
///---------------------------------------------------------------------------//
/// Material has only one component,
/// Density is either object density or it is ignored ==> use object density
/// This function is not called for a compound, but it is called on the
/// material-ized components of compounds
std::string pyne::Material::fluka_material_component(int fid, int nucid,
                                               std::string fluka_name) {
  int znum = pyne::nucname::znum(nucid);

  double atomic_mass;
  if (0 != pyne::NUC_DATA_PATH.length()) {
    // for compounds (i.e., unrecognized nucids), this will be 0
    atomic_mass = pyne::atomic_mass(nucid);
  } else {
    atomic_mass = 1.0;
  }

  return fluka_material_line(znum, atomic_mass, fid, fluka_name);
}

///---------------------------------------------------------------------------//
/// fluka_material_line
///---------------------------------------------------------------------------//
/// Given all the info, return the Material string
std::string pyne::Material::fluka_material_line(int znum, double atomic_mass,
                                          int fid, std::string fluka_name) {
  std::stringstream ls;

  if (metadata.isMember("comments") ) {
     std::string comment = metadata["comments"].asString();
     ls << "* " << comment;
     ls << std::endl;
  }
  ls << std::setw(10) << std::left << "MATERIAL";
  ls << std::setprecision(0) << std::fixed << std::showpoint <<
        std::setw(10) << std::right << (float)znum;

  ls << fluka_format_field(atomic_mass);
  // Note this is the current object density, and may or may not be meaningful
  ls << fluka_format_field(std::sqrt(density*density));

  ls << std::setprecision(0) << std::fixed << std::showpoint <<
        std::setw(10) << std::right << (float)fid;
  ls << std::setw(10) << std::right << "";
  ls << std::setw(10) << std::right << "";
  ls << std::setw(10) << std::left << fluka_name << std::endl;

  return ls.str();
}

///---------------------------------------------------------------------------//
/// fluka_format_field
///---------------------------------------------------------------------------//
/// Convenience function that returns a 10-character formatted string
/// 999 -> 999.
/// 999.12 -> 999.12
/// 999.123 -> 999.123
/// 999.1234 -> 999.123
std::string pyne::Material::fluka_format_field(float field) {
  std::stringstream ls;
  double intpart;
  modf (field, &intpart);
  if (field == intpart) {
    ls << std::setprecision(0) << std::fixed << std::showpoint
       << std::setw(10) << std::right << field;
  } else {
  // This will print however many digits after the decimal, up to a max of six
    ls.unsetf(std::ios::showpoint);
    ls.unsetf(std::ios::floatfield);
    ls.precision(6);
    ls << std::setw(10) << std::right << field;
  }

  return ls.str();
}

///---------------------------------------------------------------------------//
/// fluka_compound_str
///---------------------------------------------------------------------------//
/// Returns
/// -- MATERIAL line for compound
/// -- COMPOUND lines
std::string pyne::Material::fluka_compound_str(int id, std::string frac_type) {
  std::stringstream ss;
  std::map<double, std::string> frac_name_map;
  std::string compound_string = "";
  std::vector<std::string> material_names;

  // The nucid doesn't make sense for a compound
  int znum = 1;
  double atomic_mass = 1.;
  // This better be true
  std::string compound_name;
  if (metadata.isMember("fluka_name")) {
    compound_name = metadata["fluka_name"].asString();
  } else {
    std::cerr << "Error:  metadata \"fluka_name\" expected." << std::endl;
    compound_name = "NotFound";
  }
  ss << fluka_material_line(znum, atomic_mass, id, compound_name);

  std::string frac_sign;
  if ("atom" == frac_type) {
    frac_sign = "";
  } else {
    frac_sign = "-";
  }

  std::stringstream temp_s;
  temp_s << std::scientific;
  temp_s << std::setprecision(3);

  int counter = comp.size();
  pyne::comp_iter nuc = comp.begin();
  // This will pick up multiples of 3 components
  while (counter >= 3) {
    ss << std::setw(10) << std::left  << "COMPOUND";

    temp_s << frac_sign << nuc->second;

    ss << std::setw(10) << std::right << temp_s.str();
    ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
    nuc++;
    temp_s.str("");  // reset the stringstream for reuse

    temp_s << frac_sign << nuc->second;
    ss << std::setw(10) << std::right << temp_s.str();
    ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
    nuc++;
    temp_s.str("");

    temp_s << frac_sign << nuc->second;
    ss << std::setw(10) << std::right << temp_s.str();
    ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
    nuc++;
    temp_s.str("");

    ss << std::setw(10) << std::left << compound_name;
    ss << std::endl;

    counter -= 3;
  }

  // Get the last (or only, as the case may be) one or two fractions
  if (nuc != comp.end()) {
    ss << std::setw(10) << std::left  << "COMPOUND";
    temp_s << frac_sign << nuc->second;
    ss << std::setw(10) << std::right << temp_s.str();
    ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
    nuc++;
    temp_s.str("");

    if  (nuc != comp.end()) {
      temp_s << frac_sign << nuc->second;
      ss << std::setw(10) << std::right << temp_s.str();
      ss << std::setw(10) << std::right << nucname::fluka(nuc->first);
      nuc++;
      temp_s.str("");
    } else {
      ss << std::setw(10) << std::right << "";
      ss << std::setw(10) << std::right << "";
    }

    ss << std::setw(10) << std::right << "";
    ss << std::setw(10) << std::right << "";
    ss << std::setw(10) << std::left << compound_name;
    ss << std::endl;
    }

  return ss.str();
}

void pyne::Material::from_text(char * filename) {
  std::string fname (filename);
  from_text(fname);
}


void pyne::Material::from_text(std::string filename) {
  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // New filestream
  std::ifstream f;
  f.open(filename.c_str());

  // Read in
  comp.clear();
  std::string keystr, valstr;

  f >> keystr;
  while ( !f.eof() ) {

    if (0 == keystr.length())
      continue;

    if (keystr == "Mass"){
      f >> valstr;
      mass = pyne::to_dbl(valstr);
    } else if (keystr == "Density") {
      f >> valstr;
      density = pyne::to_dbl(valstr);
    } else if (keystr == "APerM") {
      f >> valstr;
      atoms_per_molecule = pyne::to_dbl(valstr);
    } else if (pyne::nucname::isnuclide(keystr) ||
               pyne::nucname::iselement(keystr)) {
      f >> valstr;
      if (comp.count(pyne::nucname::id(keystr))>0) {
        comp[pyne::nucname::id(keystr)] += pyne::to_dbl(valstr);
      } else {
        comp[pyne::nucname::id(keystr)] = pyne::to_dbl(valstr);
      }
    } else {
      getline(f, valstr);
      valstr= valstr.substr(0, valstr.length()-1);
      metadata[keystr]= valstr;

    }
    f >> keystr;
   }

   f.close();
   norm_comp();
}



void pyne::Material::write_text(char * filename) {
  std::string fname (filename);
  write_text(fname);
}


void pyne::Material::write_text(std::string filename) {
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);

  Json::Reader reader;
  std::vector<std::string> obj = metadata.getMemberNames();

  if (0 <= mass)
    f << "Mass    " << mass << "\n";

  if (0 <= density)
    f << "Density "  << density << "\n";

  if (0 <= atoms_per_molecule)
    f << "APerM   " << atoms_per_molecule << "\n";

  for (int i=0; i < metadata.size(); i=i+2){
    f <<metadata.get(obj.at(i), "") << metadata.get(obj.at(i+1), "");
  }

  std::string nuc_name;
  for(pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    nuc_name = pyne::nucname::name( i->first ) + "  ";
    while (nuc_name.length() < 8)
      nuc_name += " ";
    f << nuc_name << i->second << "\n";
  }

  f.close();
}


void pyne::Material::load_json(Json::Value json) {
  Json::Value::Members keys = json["comp"].getMemberNames();
  Json::Value::Members::const_iterator ikey = keys.begin();
  Json::Value::Members::const_iterator ikey_end = keys.end();
  comp.clear();
  for (; ikey != ikey_end; ++ikey)
    comp[nucname::id(*ikey)] = json["comp"][*ikey].asDouble();
  norm_comp();
  mass = json["mass"].asDouble();
  density = json["density"].asDouble();
  atoms_per_molecule = json["atoms_per_molecule"].asDouble();
  metadata = json["metadata"];
}


Json::Value pyne::Material::dump_json() {
  Json::Value json = Json::Value(Json::objectValue);
  Json::Value jcomp = Json::Value(Json::objectValue);
  json["mass"] = mass;
  json["density"] = density;
  json["atoms_per_molecule"] = atoms_per_molecule;
  json["metadata"] = metadata;
  for(comp_iter i = comp.begin(); i != comp.end(); i++)
    jcomp[nucname::name(i->first)] = (i->second);
  json["comp"] = jcomp;
  return json;
}


void pyne::Material::from_json(char * filename) {
  std::string fname (filename);
  from_json(fname);
}

void pyne::Material::from_json(std::string filename) {
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);
  std::string s;
  std::ifstream f (filename.c_str(), std::ios::in | std::ios::binary);
  f.seekg(0, std::ios::end);
  s.resize(f.tellg());
  f.seekg(0, std::ios::beg);
  f.read(&s[0], s.size());
  f.close();
  Json::Reader reader;
  Json::Value json;
  reader.parse(s, json);
  load_json(json);
}


void pyne::Material::write_json(char * filename) {
  std::string fname (filename);
  write_json(fname);
}

void pyne::Material::write_json(std::string filename) {
  Json::Value json = dump_json();
  Json::StyledWriter writer;
  std::string s = writer.write(json);
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);
  f << s << "\n";
  f.close();
}


/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/

pyne::Material::Material() {
  // Empty Material constructor
  mass = -1.0;
  density = -1.0;
  atoms_per_molecule = -1.0;
  metadata = Json::Value(Json::objectValue);
}


pyne::Material::Material(pyne::comp_map cm, double m, double d, double apm,
                         Json::Value attributes) {
  // Initializes the mass stream based on an isotopic component dictionary.
  comp = cm;
  mass = m;
  density=d;
  atoms_per_molecule = apm;
  metadata = attributes;
  if (!comp.empty())
    norm_comp();
}



pyne::Material::Material(char * filename, double m, double d, double apm,
                         Json::Value attributes) {
  mass = m;
  density=d;
  atoms_per_molecule = apm;
  metadata = attributes;

  // Check that the file is there
  std::string fname (filename);
  if (!pyne::file_exists(fname))
    throw pyne::FileNotFound(fname);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(fname.c_str());
  if (ish5)
    from_hdf5(fname);
  else
    from_text(fname);
}


pyne::Material::Material(std::string filename, double m, double d, double apm,
                         Json::Value attributes) {
  // Initializes the mass stream based on an isotopic composition file with a string name.
  mass = m;
  density=d;
  atoms_per_molecule = apm;
  metadata = attributes;

  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // Check to see if the file is in HDF5 format.
  bool ish5 = H5Fis_hdf5(filename.c_str());
  if (ish5)
    from_hdf5(filename);
  else
    from_text(filename);
}


pyne::Material::~Material() {
}



/*--- Method definitions ---*/


std::ostream& operator<<(std::ostream& os, pyne::Material mat) {
  //print the Mass Stream to stdout
  os << "\tMass: " << mat.mass << "\n";
  os << "\t---------\n";
  for(pyne::comp_iter i = mat.comp.begin(); i != mat.comp.end(); i++)
  {
    os << "\t" << pyne::nucname::name( i->first ) << "\t" << i->second << "\n";
  }
  return os;
}

// Note this refines << for an inheritor of std::ostream.
std::ostringstream& operator<<(std::ostringstream& os, pyne::Material mat) {
  return os;
}

void pyne::Material::normalize () {
  // normalizes the mass
  mass = 1.0;
}


pyne::comp_map pyne::Material::mult_by_mass() {
  // bypass calculation if already normalized.
  if (mass == 1.0)
    return comp;

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    cm[i->first] = (i->second) * mass;
  }
  return cm;
}



pyne::comp_map pyne::Material::activity() {
  pyne::comp_map act;
  double masspermole = mass * pyne::N_A;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
    act[i->first] = masspermole * (i->second) * decay_const(i->first) / \
                    atomic_mass(i->first);
  }
  return act;
}


pyne::comp_map pyne::Material::decay_heat() {
  pyne::comp_map dh;
  double masspermole = mass * pyne::N_A;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
    dh[i->first] = masspermole * (i->second) * \
                   decay_const(metastable_id(i->first,nucname::snum(i->first))) * \
                   q_val(i->first) / atomic_mass(i->first) / pyne::MeV_per_MJ;
  }
  return dh;
}


pyne::comp_map pyne::Material::dose_per_g(std::string dose_type, int source) {
  pyne::comp_map dose;
  const double pCi_per_Bq = 27.027027;
  if (dose_type == "ext_air") {
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
      dose[i->first] = Ci_per_Bq * pyne::N_A * (i->second) * \
                       decay_const(i->first) * ext_air_dose(i->first, source) / \
                       atomic_mass(i->first);
    }
  } else if (dose_type == "ext_soil") {
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
      dose[i->first] = Ci_per_Bq * pyne::N_A * (i->second) * \
                       decay_const(i->first) * ext_soil_dose(i->first, source) / \
                       atomic_mass(i->first);
    }
  } else if (dose_type == "ingest") {
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
      dose[i->first] = pCi_per_Bq * pyne::N_A * (i->second) * \
                       decay_const(i->first) * ingest_dose(i->first, source) / \
                       atomic_mass(i->first);
    }
  } else if (dose_type == "inhale") {
    for (pyne::comp_iter i = comp.begin(); i != comp.end(); ++i) {
      dose[i->first] = pCi_per_Bq * pyne::N_A * (i->second) * \
                       decay_const(i->first) * inhale_dose(i->first, source) / \
                       atomic_mass(i->first);
    }
  } else {
    throw std::invalid_argument("Dose type must be one of: ext_air, ext_soil, ingest, inhale.");
  }
  return dose;
}


double pyne::Material::molecular_mass(double apm) {
  // Calculate the atomic weight of the Material
  double inverseA = 0.0;

  for (pyne::comp_iter nuc = comp.begin(); nuc != comp.end(); nuc++)
    inverseA += (nuc->second) / pyne::atomic_mass(nuc->first);

  if (inverseA == 0.0)
    return inverseA;

  // select the atoms per mol
  double atsperm = 1.0; // default to 1.0
  if (0.0 <= apm) {
    atsperm = apm;            // take the function argument, if valid
    if (atoms_per_molecule < 0.0)
      atoms_per_molecule = apm;     // Store the function argument on class, if class has no value
  } else if (0.0 <= atoms_per_molecule)
    atsperm = atoms_per_molecule;  // select the class's value

  return atsperm / inverseA;
}

pyne::Material pyne::Material::expand_elements(std::set<int> exception_ids) {
  // Expands the natural elements of a material and returns a new material note
  // that this implementation relies on the fact that maps of ints are stored in
  // a sorted manner in C++.
  int n, nabund, znuc, zabund;
  comp_map newcomp;
  std::map<int, double>::iterator abund_itr, abund_end;
  if (pyne::natural_abund_map.empty())
    pyne::_load_atomic_mass_map();
  abund_itr = pyne::natural_abund_map.begin();
  abund_end = pyne::natural_abund_map.end();
  zabund = nucname::znum((*abund_itr).first);
  for (comp_iter nuc = comp.begin(); nuc != comp.end(); nuc++) {
    // keep element as-is if in exception list
    if (0 < exception_ids.count(nuc->first)) {
      newcomp.insert(*nuc);
      continue;
    }

    if(abund_itr == abund_end)
      newcomp.insert(*nuc);
    else if(0 == nucname::anum((*nuc).first)) {
      n = (*nuc).first;
      znuc = nucname::znum(n);
      if (znuc < zabund) {
        newcomp.insert(*nuc);
        continue;
      }
      while(zabund <= znuc) {
        nabund = (*abund_itr).first;
        zabund = nucname::znum(nabund);

        if (zabund == znuc && 0 != nucname::anum(nabund) && 0.0 != (*abund_itr).second)
          newcomp[nabund] = (*abund_itr).second * (*nuc).second * \
                            atomic_mass(nabund) / atomic_mass(n);
        else if (n == nabund && 0.0 == (*abund_itr).second)
          newcomp.insert(*nuc);
        abund_itr++;
        if (abund_itr == abund_end) {
          zabund = INT_MAX;
          break;
        }
      }
    } else
      newcomp.insert(*nuc);
  }
  return Material(newcomp, mass, density, atoms_per_molecule, metadata);
}

// Wrapped version for calling from python
pyne::Material pyne::Material::expand_elements(int** int_ptr_arry ) {
    std::set<int> nucvec;
    // Set first pointer to first int pointed to by arg
    if (int_ptr_arry != NULL) {
      int *int_ptr = *int_ptr_arry;
      while (int_ptr != NULL)
  {
    nucvec.insert(*int_ptr);
    int_ptr++;
  }
    }
    return expand_elements(nucvec);
}


pyne::Material pyne::Material::collapse_elements(std::set<int> exception_ids) {
  ////////////////////////////////////////////////////////////////////////
  // Assumptions
  //    - list passed in is of nucid's formed from the znum-anum of
  //      Fluka-named isotopes, since we want to preserve the full
  //      nucid of any such material in the problem
  // Algorithm
  // for each component listed in this material that has a nonzero frac or
  //    weight amount, look at its 'stripped' nucid, that is, the last four
  //    places replaced by zeros.
  //    if it's on the exception list, copy the component
  //    else it is to be collapsed
  //       => add its frac to the component of the znum
  //
  // * When from_hdf5 reads from a file the comp iterator will produce a
  //   hit for EVERY nucid in EVERY material in the file.  Only the nucids
  //   belonging to the CURRENT material have a nonzero fraction/mass amount
  /////////////////////////////////////////////////////////////////////////
  pyne::comp_map cm;

  for (pyne::comp_iter ptr = comp.begin(); ptr != comp.end(); ptr++) {
      if (0 < ptr->second) {
        // There is a nonzero amount of this nucid in the current material,
        // check if znum and anum are in the exception list,
        int cur_stripped_id = nucname::znum(ptr->first)*10000000
                        + nucname::anum(ptr->first)*10000;
        if (0 < exception_ids.count(cur_stripped_id)) {
        // The znum/anum combination identify the current material as a
        // fluka-named exception list => copy, don't collapse
          cm[ptr->first] = (ptr->second) * mass;
        } else {
          // Not on exception list => add frac to id-component
          int znum_id = nucname::id(nucname::znum(ptr->first));
          cm[znum_id] += (ptr->second) * mass;
        }
      }
  }
  // Copy
  pyne::Material collapsed = pyne::Material(cm, mass, density,
                                            atoms_per_molecule, metadata);
  return collapsed;
}

// Wrapped version for calling from python
pyne::Material pyne::Material::collapse_elements(int** int_ptr_arry ) {
    std::set<int> nucvec;
    // Set first pointer to first int pointed to by arg
    int *int_ptr = *int_ptr_arry;
    while (int_ptr != NULL)
    {
      nucvec.insert(*int_ptr);
      int_ptr++;
    }
    return collapse_elements(nucvec);
}

  // Set up atom or mass frac map

std::map<int, double> pyne::Material::get_density_frac(std::string frac_type){
  std::map<int, double> fracs;

  if ("atom" == frac_type) {
    if (density != -1.0) {
      fracs = to_atom_dens();
      for (comp_iter ci = fracs.begin(); ci != fracs.end(); ci++){
        ci->second *= pyne::cm2_per_barn; // unit requirememt is [10^24 atoms/cm3] = [atoms/b.cm]
      }
    } else {
      fracs = to_atom_frac();
    }
  } else {
    fracs = comp;
    if (density != -1.0) {
      for (comp_iter ci = fracs.begin(); ci != fracs.end(); ci++){
        ci->second *= density;
      }
    }
  }
  return fracs;
}



double pyne::Material::mass_density(double num_dens, double apm) {
  if (0.0 <= num_dens) {
    double mw = molecular_mass(apm);
    density = num_dens * mw / pyne::N_A / atoms_per_molecule;
  }
  return density;
}


double pyne::Material::number_density(double mass_dens, double apm) {
  if (0 <= mass_dens)
    density = mass_dens;
  double mw = molecular_mass(apm);
  double num_dens = density * pyne::N_A * atoms_per_molecule / mw;
  return num_dens;
}


/*--- Stub-Stream Computation ---*/

pyne::Material pyne::Material::sub_mat(std::set<int> nucset) {
  // Grabs a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ( 0 < nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  }

  return pyne::Material(cm, -1, -1);
}



pyne::Material pyne::Material::sub_mat(std::set<std::string> nucset) {
  // Grabs a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  }

  return sub_mat(iset);
}



pyne::Material pyne::Material::set_mat (std::set<int> nucset, double value) {
  // Sets a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;

  // Add non-set components
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ( 0 == nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  }

  // Add set component
  for (std::set<int>::iterator nuc = nucset.begin(); nuc != nucset.end(); nuc++)
    cm[*nuc] = value;

  return pyne::Material(cm, -1, -1);
}



pyne::Material pyne::Material::set_mat(std::set<std::string> nucset, double value) {
  // Sets a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  }

  return set_mat(iset, value);
}




pyne::Material pyne::Material::del_mat(std::set<int> nucset) {
  // Removes a sub-material from this mat based on a set of integers.
  // Integers can either be of id form -OR- they can be a z-numer (is 8 for O, 93 for Np, etc).
  // n is the name of the new material.

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    // Only add to new comp if not in nucset
    if ( 0 == nucset.count(i->first) )
      cm[i->first] = (i->second) * mass;
  }

  return pyne::Material(cm, -1, -1);
}



pyne::Material pyne::Material::del_mat (std::set<std::string> nucset) {
  // Removes a substream from this stream based on a set of strings.
  // Strings can be of any form.
  std::set<int> iset;
  for (std::set<std::string>::iterator i = nucset.begin(); i != nucset.end(); i++) {
    iset.insert(pyne::nucname::id(*i));
  }

  return del_mat(iset);
}






pyne::Material pyne::Material::sub_range(int lower, int upper) {
  // Grabs a sub-material from this mat based on a range of integers.
  if (upper < lower)
  {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  }

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ((lower <= (i->first)) && ((i->first) < upper))
      cm[i->first] = (i->second) * mass;
  }

  return pyne::Material(cm, -1,-1);
}



pyne::Material pyne::Material::set_range(int lower, int upper, double value) {
  // Sets a sub-material from this mat based on a range of integers.
  if (upper < lower) {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  }

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ((lower <= (i->first)) && ((i->first) < upper))
      cm[i->first] = value;
    else
      cm[i->first] = (i->second) * mass;
  }

  return pyne::Material(cm, -1,-1);
}



pyne::Material pyne::Material::del_range(int lower, int upper) {
  // Removes a sub-material from this mat based on a range of integers.
  if (upper < lower) {
    int temp_upper = upper;
    upper = lower;
    lower = temp_upper;
  }

  pyne::comp_map cm;
  for (pyne::comp_iter i = comp.begin(); i != comp.end(); i++) {
    if ((upper <= (i->first)) || ((i->first) < lower))
      cm[i->first] = (i->second) * mass;
  }

  return pyne::Material(cm, -1, -1);
}










pyne::Material pyne::Material::sub_elem(int elem) {
  // Returns a material of the element that is a submaterial of this one.
  return sub_range(elem, elem + 10000000);
}



pyne::Material pyne::Material::sub_lan() {
  // Returns a material of Lanthanides that is a sub-material of this one.
  return sub_range(570000000, 720000000);
}



pyne::Material pyne::Material::sub_act() {
  //Returns a material of Actindes that is a sub-material of this one.
  return sub_range(890000000, 1040000000);
}


pyne::Material pyne::Material::sub_tru() {
  // Returns a material of Transuranics that is a sub-material of this one.
  return sub_range(930000000, INT_MAX);
}



pyne::Material pyne::Material::sub_ma() {
  // Returns a material of Minor Actinides that is a sub-material of this one.
  return sub_range(930000000, 1040000000).del_range(940000000, 950000000);
}



pyne::Material pyne::Material::sub_fp() {
  // Returns a material of Fission Products that is a sub-material of this one.
  return sub_range(0, 890000000);
}




/*--- Atom Frac Functions ---*/

std::map<int, double> pyne::Material::to_atom_frac() {
  // Returns an atom fraction map from this material's composition
  // the material's molecular mass
  double mat_mw = molecular_mass();

  std::map<int, double> atom_fracs = std::map<int, double>();

  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++)
    atom_fracs[ci->first] = (ci->second) * mat_mw / pyne::atomic_mass(ci->first);

  return atom_fracs;
}


void pyne::Material::from_atom_frac(std::map<int, double> atom_fracs) {
  // atom frac must be of the form {nuc: af}, eg, water
  //  80160: 1.0
  //  10010: 2.0

  // clear existing components
  comp.clear();
  atoms_per_molecule = 0.0;

  for (std::map<int, double>::iterator afi = atom_fracs.begin(); afi != atom_fracs.end(); afi++) {
    comp[afi->first] = (afi->second) * pyne::atomic_mass(afi->first);
    atoms_per_molecule += (afi->second);
  }

  norm_comp();
}


std::map<int, double> pyne::Material::to_atom_dens() {
  // Returns an atom density map from this material's composition
  // the material's density

  std::map<int, double> atom_dens = std::map<int, double>();

  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++)
    atom_dens[ci->first] = (ci->second) * density * pyne::N_A / pyne::atomic_mass(ci->first);

  return atom_dens;
}


std::vector<std::pair<double, double> > pyne::Material::gammas() {
  std::vector<std::pair<double, double> > result;
  std::map<int, double> atom_fracs = this->to_atom_frac();
  int state_id;
  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++) {
    if (ci->first % 10000 > 0)
        state_id = nucname::id_to_state_id(ci->first);
    else
        state_id = ci->first;

    std::vector<std::pair<double, double> > raw_gammas = pyne::gammas(state_id);
    for (int i = 0; i < raw_gammas.size(); ++i) {
      result.push_back(std::make_pair(raw_gammas[i].first,
        atom_fracs[ci->first]*raw_gammas[i].second));
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::Material::xrays() {
  std::vector<std::pair<double, double> > result;
  std::map<int, double> atom_fracs = this->to_atom_frac();
  int state_id;
  for (comp_iter ci = comp.begin(); ci != comp.end(); ci++) {
    if (ci->first % 10000 > 0)
        state_id = nucname::id_to_state_id(ci->first);
    else
        state_id = ci->first;

    std::vector<std::pair<double, double> > raw_xrays = pyne::xrays(state_id);
    for (int i = 0; i < raw_xrays.size(); ++i) {
      result.push_back(std::make_pair(raw_xrays[i].first,
        atom_fracs[ci->first]*raw_xrays[i].second));
    }
  }
  return result;
}

std::vector<std::pair<double, double> > pyne::Material::photons(bool norm) {
  std::vector<std::pair<double, double> >  txray = this->xrays();
  std::vector<std::pair<double, double> >  tgammas = this->gammas();
  for (int i = 0; i < txray.size(); ++i)
    tgammas.push_back(txray[i]);
  if (norm)
    tgammas = normalize_radioactivity(tgammas);
  return tgammas;
}

std::vector<std::pair<double, double> > pyne::Material::normalize_radioactivity(
std::vector<std::pair<double, double> > unnormed) {
  std::vector<std::pair<double, double> > normed;
  double sum = 0.0;
  for (int i = 0; i < unnormed.size(); ++i) {
    if (!isnan(unnormed[i].second))
      sum = sum + unnormed[i].second;
  }
  for (int i = 0; i < unnormed.size(); ++i) {
    if (!isnan(unnormed[i].second)) {
      normed.push_back(std::make_pair(unnormed[i].first,
        (unnormed[i].second)/sum));
    }
  }
  return normed;
}


pyne::Material pyne::Material::decay(double t) {
  Material rtn;
  comp_map out = pyne::decayers::decay(to_atom_frac(), t);
  rtn.from_atom_frac(out);
  rtn.mass = mass * rtn.molecular_mass() / molecular_mass();
  return rtn;
}


pyne::Material pyne::Material::cram(std::vector<double> A,
                                    const int order) {
  Material rtn;
  rtn.from_atom_frac(pyne::transmuters::cram(A, to_atom_frac(), order));
  rtn.mass = mass * rtn.molecular_mass() / molecular_mass();
  return rtn;
}

pyne::Material pyne::Material::operator+ (double y) {
  // Overloads x + y
  return pyne::Material(comp, mass + y, density);
}



pyne::Material pyne::Material::operator+ (Material y) {
  // Overloads x + y
  pyne::comp_map cm;
  pyne::comp_map xwgt = mult_by_mass();
  pyne::comp_map ywgt = y.mult_by_mass();

  for (pyne::comp_iter i = xwgt.begin(); i != xwgt.end(); i++) {
    if ( 0 < ywgt.count(i->first) )
      cm[i->first] = xwgt[i->first] + ywgt[i->first];
    else
      cm[i->first] = xwgt[i->first];
  }

  for (pyne::comp_iter i = ywgt.begin(); i != ywgt.end(); i++) {
    if ( 0 == cm.count(i->first) )
      cm[i->first] = ywgt[i->first];
  }

  return pyne::Material(cm, -1, -1);
}



pyne::Material pyne::Material::operator* (double y) {
  // Overloads x * y
  return pyne::Material(comp, mass * y, density);
}



pyne::Material pyne::Material::operator/ (double y) {
  // Overloads x / y
  return pyne::Material(comp, mass / y, density );
}
//
// end of /home/mouginot/work/app/pyne/src/material.cpp
//


//
// start of /home/mouginot/work/app/pyne/src/material_library.cpp
//
#include <unistd.h>
#include <iostream>

#ifndef PYNE_IS_AMALGAMATED
#include "material_library.h"
#endif

// Empty Constructor
pyne::MaterialLibrary::MaterialLibrary(){};

// Default constructor
pyne::MaterialLibrary::MaterialLibrary(const std::string& file,
                                       const std::string& datapath,
                                       const std::string& nucpath) {
  if (!pyne::file_exists(file)) {
    throw std::runtime_error("File " + file +
                             " not found or no read permission");
  }
  if (!hdf5_path_exists(file, datapath)) {
    throw std::runtime_error("The datapath, " + datapath + ", in " + file +
                             " is empty.");
  }
  // load materials
  from_hdf5(file, datapath, nucpath);
};

void pyne::MaterialLibrary::from_hdf5(char* fname, char* dpath, char* npath,
                                      int protocol) {
  std::string filename(fname);
  std::string datapath(dpath);
  std::string nucpath(npath);
  from_hdf5(filename, datapath, npath, protocol);
}

// Append Material to the library from hdf5 file
void pyne::MaterialLibrary::from_hdf5(const std::string& filename,
                                      const std::string& datapath,
                                      const std::string& nucpath,
                                      int protocol) {
  if (!hdf5_path_exists(filename, datapath)){
    throw std::runtime_error("The datapath, " + datapath + ", in " + filename +
                             " is empty.");
  }

  int file_num_materials = get_length_of_table(filename, datapath);
  int library_length = material_library.size();
  for (int i = 0; i < file_num_materials; i++) {
    pyne::Material mat = pyne::Material();  
    // Get material from the hdf5 file
    if (nucpath == "") {
      mat.from_hdf5(filename, datapath, i, protocol);
    } else {
      mat.from_hdf5(filename, datapath, nucpath, i, protocol);
    }
    (*this).add_material(mat);
  }
}

void pyne::MaterialLibrary::merge(pyne::MaterialLibrary mat_lib) {
  pyne::matname_set mats_to_add = mat_lib.get_keylist();
  for (auto it = mats_to_add.begin(); it != mats_to_add.end(); it++) {
    pyne::Material mat = Material(mat_lib.get_material(*it));
    (*this).add_material(*it, mat);
  }
}

void pyne::MaterialLibrary::load_json(Json::Value json) {
  Json::Value::Members keys = json.getMemberNames();
  Json::Value::Members::const_iterator ikey = keys.begin();
  Json::Value::Members::const_iterator ikey_end = keys.end();
  for (; ikey != ikey_end; ++ikey) {
    pyne::Material mat = pyne::Material();
    mat.load_json(json[*ikey]);
    (*this).add_material(*ikey, Material(mat));
  }
}

Json::Value pyne::MaterialLibrary::dump_json() {
  Json::Value json = Json::Value(Json::objectValue);

  for (auto name : name_order) {
    json[name] = material_library[name]->dump_json();
  }
  return json;
}

void pyne::MaterialLibrary::from_json(char* fname) {
  std::string filename(fname);
  from_json(filename);
}
void pyne::MaterialLibrary::from_json(const std::string& filename) {
  if (!pyne::file_exists(filename)) throw pyne::FileNotFound(filename);
  std::string s;
  std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
  f.seekg(0, std::ios::end);
  s.resize(f.tellg());
  f.seekg(0, std::ios::beg);
  f.read(&s[0], s.size());
  f.close();
  Json::Reader reader;
  Json::Value json;
  reader.parse(s, json);
  load_json(json);
}

void pyne::MaterialLibrary::write_json(char* filename) {
  std::string fname(filename);
  write_json(fname);
}

void pyne::MaterialLibrary::write_json(const std::string& filename) {
  Json::Value json = dump_json();
  Json::StyledWriter writer;
  std::string s = writer.write(json);
  std::ofstream f;
  f.open(filename.c_str(), std::ios_base::trunc);
  f << s << "\n";
  f.close();
}

void pyne::MaterialLibrary::add_material(pyne::Material mat) {
  // if exists, get the material name from metadata make one instead
  std::string mat_name;
  int mat_number = -1;
  std::set<int>::iterator mat_numb_it;
  if (mat.metadata.isMember("mat_number")) {
    if (Json::intValue <= mat.metadata["mat_number"].type() &&
        mat.metadata["mat_number"].type() <= Json::realValue ){
      mat_number = mat.metadata["mat_number"].asInt();
    } else {
      mat_number = std::stoi(mat.metadata["mat_number"].asString());
    }
  } 
  pyne::matname_set::iterator key_it;
  if (mat.metadata.isMember("name")) {
    
    if (Json::intValue <= mat.metadata["name"].type() &&
        mat.metadata["name"].type() <= Json::realValue ){
      mat_name = std::to_string(mat.metadata["name"].asInt());
      mat.metadata["name"] = mat_name; 
    } else {
      mat_name = mat.metadata["name"].asString();
    }
    
  } else {
    if (mat_number == -1) {
      mat_number = keylist.size(); //set a temp mat_number to form mat name 
    }
    mat_name = "_" + std::to_string(mat_number);
    mat.metadata["name"] = mat_name;
  }
  
  add_material(mat_name, mat);  
}

void pyne::MaterialLibrary::add_material(char* key, pyne::Material mat) {
  std::string key_name(key);
  add_material(key_name, mat);
}

void pyne::MaterialLibrary::add_material(const std::string& key, pyne::Material mat) {
  std::string mat_name;
  int mat_number = -1;
  std::set<int>::iterator mat_numb_it;
  if (mat.metadata.isMember("mat_number")) {
    if (Json::intValue <= mat.metadata["mat_number"].type() &&
        mat.metadata["mat_number"].type() <= Json::realValue ){
      mat_number = mat.metadata["mat_number"].asInt();
    } else {
      mat_number = std::stoi(mat.metadata["mat_number"].asString());
      mat.metadata["mat_number"] = mat_number;
    }
    mat_numb_it = mat_number_set.find(mat_number);
    if (mat_numb_it != mat_number_set.end()) {
      std::string msg = "Material number ";
      msg += mat_number;
      msg += " is already in the library.";
      warning(msg);
    }
  }  
  if ( !mat.metadata.isMember("name")) {
    mat.metadata["name"] = key;
  } 

  append_to_nuclist(mat);
  if(mat_number > 0)
    mat_number_set.insert(mat_number);
  std::pair<std::set<std::string>::iterator, bool> key_insert;
  key_insert = keylist.insert(key);
  // if the key was not present in the MaterialLibrary add it, otherwise overwrite the existing material
  if( key_insert.second == true)
    name_order.push_back(key);
  material_library[key] = new Material(mat);
}

void pyne::MaterialLibrary::del_material(pyne::Material mat) {
  if (mat.metadata.isMember("name")) {
    std::string mat_name = mat.metadata["name"].asString();
    del_material(mat_name);
  }
}

void pyne::MaterialLibrary::del_material(const std::string& key) {
  material_library.erase(key);
  keylist.erase(key);
  for (auto name = name_order.begin(); name != name_order.end(); name++) {
    if (*name == key) {
      name_order.erase(name);
      return;
    }
  }
}

pyne::Material pyne::MaterialLibrary::get_material(
    const std::string& mat_name) const {
  auto it = material_library.find(mat_name);
  if (it != material_library.end()) {
    return *(it->second);
  } else {
    return pyne::Material();
  }
}

pyne::Material* pyne::MaterialLibrary::get_material_ptr(
    const std::string& mat_name) {
  auto it = material_library.find(mat_name);
  if (it != material_library.end()) {
    return it->second;
  } else {
    return new pyne::Material();
  }
}

void pyne::MaterialLibrary::replace(int num, pyne::Material mat) {
  if (!mat.metadata.isMember("name")) {
    mat.metadata["name"] = name_order[num];
  } else if (mat.metadata["name"] != name_order[num]) {
    material_library.erase(name_order[num]);
    name_order[num] = mat.metadata["name"].asString();
  }
  material_library[name_order[num]] = new pyne::Material(mat);
}

void pyne::MaterialLibrary::write_hdf5(char* fname, char* dpath, char* npath) {
  std::string filename(fname);
  std::string datapath(dpath);
  std::string nucpath(npath);
  write_hdf5(filename, datapath, nucpath);
}

void pyne::MaterialLibrary::write_hdf5(const std::string& filename,
                                       const std::string& datapath,
                                       const std::string& nucpath) {
  // A large part of this is inspired/taken from by material.cpp...
  // Turn off annoying HDF5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl;
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  // Create new/open datafile.
  hid_t db;
  if (pyne::file_exists(filename)) {
    bool ish5 = H5Fis_hdf5(filename.c_str());
    if (!ish5) throw h5wrap::FileNotHDF5(filename);
    db = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
  } else
    db = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

  std::vector<int> nuclides;
  std::copy(nuclist.begin(), nuclist.end(),
            inserter(nuclides, nuclides.begin()));

  int nuclide_size;
  hsize_t nuclide_dims[1];
  nuclide_size = nuclides.size();

  // Create the nuc data
  int nuclide_data[nuclide_size];
  int n = 0;
  for (auto nuclide : nuclides) {
    nuclide_data[n] = (nuclide);
    n++;
  }
  nuclide_dims[0] = nuclide_size;

  // Write Nuclist in the hdf5 file
  hid_t nuclide_space = H5Screate_simple(1, nuclide_dims, NULL);
  hid_t nuclide_set =
      H5Dcreate2(db, nucpath.c_str(), H5T_NATIVE_INT, nuclide_space,
                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(nuclide_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           nuclide_data);
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  // Close out the Dataset
  H5Fclose(db);

  // Write the Materials in the file
  for (auto name : name_order) {
    material_library[name]->write_hdf5(filename, datapath, nucpath);
  }
}

void pyne::MaterialLibrary::append_to_nuclist(pyne::Material mat) {
  pyne::comp_map mat_comp = mat.comp;
  for (auto nuclide : mat_comp) {
    nuclist.insert(nuclide.first);
  }
}

// see if path exists before we go on
bool pyne::MaterialLibrary::hdf5_path_exists(const std::string& filename,
                                             const std::string& datapath) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  bool datapath_exists = h5wrap::path_exists(db, datapath.c_str());
  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);

  return datapath_exists;
}

int pyne::MaterialLibrary::get_length_of_table(const std::string& filename,
                                               const std::string& datapath) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  hid_t ds = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(ds);

  hsize_t arr_dims[1];
  int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);

  return arr_dims[0];
}
//
// end of /home/mouginot/work/app/pyne/src/material_library.cpp
//


//
// start of /home/mouginot/work/app/pyne/src/particle.cpp
//
#ifndef PYNE_IS_AMALGAMATED
#include "particle.h"
#endif

std::string pyne::particle::_names[NUM_PARTICLES] = {
  // leptons
  "Electron",
  "Positron",
  "ElectronNeutrino",
  "ElectronAntiNeutrino",
  "Muon",
  "AntiMuon",
  "MuonNeutrino",
  "MuonAntiNeutrino",
  "Tauon",
  "AntiTauon",
  "TauNeutrino",
  "TauAntiNeutrino",
  // gauge bosons
  "Photon",
  // light mesons
  "Pion",
  "AntiPion",
  // strange mesons
  "Kaon",
  "AntiKaon",
  "KaonZeroShort",
  "KaonZero",
  "AntiKaonZero",
  // light baryons
  "Neutron",
  "AntiNeutron",
  "Proton",
  "AntiProton",
  // strange baryons
  "Lambda",
  "AntiLambda",
  "Sigma-",
  "AntiSigma-",
  "Sigma+",
  "AntiSigma+",
  "Sigma",
  "AntiSigmaZero"
  // Charmed baryons
};

int pyne::particle::_pdcids[NUM_PARTICLES] = {
  11,
  -11,
  12,
  -12,
  13,
  -13,
  14,
  -14,
  15,
  -15,
  16,
  -16,
  // gauge bosons
  22,
  // light mesons
  211,
  -211,
  // strange mesons
  321,
  -321,
  310,
  311,
  -311,
  // light baryons
  2112,
  -2112,
  2212,
  -2212,
  // strange Baryons
  3122,
  -3122,
  3112,
  3112,
  3222,
  -3222,
  3212,
  -3212
  // charmed baryons
};

std::set<std::string> pyne::particle::names(pyne::particle::_names,
              pyne::particle::_names+NUM_PARTICLES);

std::set<int> pyne::particle::pdc_nums(pyne::particle::_pdcids,
                 pyne::particle::_pdcids+NUM_PARTICLES);

std::map<std::string,int> pyne::particle::altnames;
std::map<int,std::string> pyne::particle::id_name;
std::map<std::string,int> pyne::particle::name_id;
std::map<std::string,std::string> pyne::particle::docs;

std::map<std::string,std::string> pyne::particle::part_to_fluka;
std::map<std::string,std::string> pyne::particle::part_to_mcnp;
std::map<std::string,std::string> pyne::particle::part_to_mcnp6;
std::map<std::string,std::string> pyne::particle::part_to_geant4;


void * pyne::particle::_fill_maps() {
  using std::make_pair;

  std::string _docs[NUM_PARTICLES] = {
    // leptons
    "Electron",
    "Positron",
    "Electron Neutrino",
    "Electron Anti Neutrino",
    "Muon Neutrino",
    "Anti Muon",
    "Muon Neutrino",
    "Muon Anti Neutrino",
    "Tauon",
    "Anti Tauon",
    "Tau Neutrino",
    "Tau Anti Neutrino",
    // gauge bosons
    "Photon",
    // light mesons
    "Pion",
    "Anti Pion",
    // strange mesons
    "Kaon",
    "Anti Kaon",
    "Kaon Zero Short",
    "Kaon Zero",
    "Anti Kaon Zero",
    // light baryons
    "Neutron",
    "Anti Neutron",
    "Proton",
    "Anti Proton",
    // strange baryons
    "Lambda",
    "Anti Lambda",
    "Sigma-",
    "Anti Sigma-",
    "Sigma+",
    "Anti Sigma+",
    "Sigma",
    "Anti Sigma Zero"
    // Charmed baryons
  };

  int pid;  // particle id
  for ( int i = 0 ; i < NUM_PARTICLES ; i++ ) {
    pid = _pdcids[i];
    // make id to name map
    id_name[pid] = _names[i];
    // make name to id map
    name_id[_names[i]] = pid;
    // make doc correspondence
    docs[_names[i]] = _docs[i];
  }

  // make the alternates
  altnames["Hydrogen"] = name_id["Proton"];
  altnames["Protium"] = name_id["Proton"];
  altnames["Beta"] = name_id["Electron"];
  altnames["Beta-"] = name_id["Electron"];
  altnames["Beta+"] = name_id["Positron"];
  altnames["Gamma"] = name_id["Photon"];
  altnames["X-Ray"] = name_id["Photon"];

  part_to_mcnp["Neutron"]="n";
  part_to_mcnp["Photon"]="p";
  part_to_mcnp["Electron"]="e";

  part_to_mcnp6["Neutron"]="n";
  part_to_mcnp6["Photon"]="p";
  part_to_mcnp6["Electron"]="e";
  part_to_mcnp6["Proton"]="h";

  part_to_fluka["Electron"]="ELECTRON";
  part_to_fluka["Positron"]="POSITRON";
  part_to_fluka["ElectronNeutrino"] ="NEUTRIE";
  part_to_fluka["ElectronAntiNeutrino"] ="ANEUTRIE";
  part_to_fluka["Muon"]="MUON+";
  part_to_fluka["AntiMuon"]="MUON-";
  part_to_fluka["MuonNeutrino"]="NEUTRIM";
  part_to_fluka["MuonAntiNeutrino"]="ANEUTRIM";
  part_to_fluka["Tauon"]="TAU+";
  part_to_fluka["Anti Tauon"]="TAU-";
  part_to_fluka["TauNeutrino"]="NEUTRIT";
  part_to_fluka["TauAntiNeutrino"]="ANEUTRIT";
  // gauge bosons
  part_to_fluka["Photon"]="PHOTON";
  // light mesons
  part_to_fluka["Pion"]="PION-";
  part_to_fluka["Anti Pion"]="PION+";
  // strange mesons
  part_to_fluka["Kaon"]="KAON+";
  part_to_fluka["AntiKaon"]="KAON-";
  part_to_fluka["KaonZero Short"]="KAONSHRT";
  part_to_fluka["KaonZero"]="KAONZERO";
  part_to_fluka["AntiKaonZero"]="AKAONZER";
  // light baryons
  part_to_fluka["Neutron"]="NEUTRON";
  part_to_fluka["AntiNeutron"]="ANEUTRON";
  part_to_fluka["Proton"]="PROTON";
  part_to_fluka["AntiProton"]="APROTON";
  // strange baryons
  part_to_fluka["Lambda"]="LAMBDA";
  part_to_fluka["AntiLambda"]="ALAMBDA";
  part_to_fluka["Sigma-"]="SIGMA-";
  part_to_fluka["Anti Sigma-"]="ASIGMA-";
  part_to_fluka["Sigma+"]="SIGMA+";
  part_to_fluka["Anti Sigma+"]="ASIGMA+";
  part_to_fluka["Sigma"]="SIGMAZER";
  part_to_fluka["AntiSigmaZero"]="ASIGMAZE";

  part_to_geant4["Electron"]="e-";
  part_to_geant4["Positron"]="e+";
  part_to_geant4["ElectronNeutrino"] ="nu_e";
  part_to_geant4["ElectronAntiNeutrino"] ="anti_nu_e";
  part_to_geant4["Muon"]="mu+";
  part_to_geant4["AntiMuon"]="mu-";
  part_to_geant4["MuonNeutrino"]="nu_mu";
  part_to_geant4["MuonAntiNeutrino"]="anti_nu_mu";
  part_to_geant4["Tauon"]="tau+";
  part_to_geant4["Anti Tauon"]="tau-";
  part_to_geant4["TauNeutrino"]="nu_tau";
  part_to_geant4["TauAntiNeutrino"]="anti_nu_tau";
  // gauge bosons
  part_to_geant4["Photon"]="gamma";
  // light mesons
  part_to_geant4["Pion"]="pi-";
  part_to_geant4["Anti Pion"]="pi+";
  // strange mesons
  part_to_geant4["Kaon"]="kaon+";
  part_to_geant4["AntiKaon"]="kaon-";
  part_to_geant4["KaonZero Short"]="kaon0S";
  part_to_geant4["KaonZero"]="kaon0";
  // light baryons
  part_to_geant4["Neutron"]="neutron";
  part_to_geant4["AntiNeutron"]="anti_neutron";
  part_to_geant4["Proton"]="proton";
  part_to_geant4["AntiProton"]="anti_proton";
  // strange baryons
  part_to_geant4["Lambda"]="lambda";
  part_to_geant4["AntiLambda"]="anti_lambda";
  part_to_geant4["Sigma-"]="sigma-";
  part_to_geant4["Anti Sigma-"]="anti_sigma-";
  part_to_geant4["Sigma+"]="sigma+";
  part_to_geant4["Anti Sigma+"]="anti_sigma+";
  part_to_geant4["Sigma"]="sigma0";
  part_to_geant4["AntiSigmaZero"]="anti_sigma0";
  return NULL;
}

void * pyne::particle::filler = pyne::particle::_fill_maps();

// is hydrogen
bool pyne::particle::is_hydrogen(int s) {
  if(s == name_id["Proton"])
    return true;
  if(pyne::particle::is_hydrogen(pyne::nucname::name(s)))
      return true;
  return false;
}

bool pyne::particle::is_hydrogen(char *s) {
  return pyne::particle::is_hydrogen(std::string(s));
}

bool pyne::particle::is_hydrogen(std::string s) {
  // check std name
  if(name_id[s] == name_id["Proton"])
    return true;
  if(altnames[s] == name_id["Proton"])
    return true;
  if(pyne::nucname::name(s).find("H1") != std::string::npos)
    return true;
  return false;
}
// heavy ion
bool pyne::particle::is_heavy_ion(int s) {
  return pyne::particle::is_heavy_ion(std::string(id_name[s]));
}

bool pyne::particle::is_heavy_ion(char *s) {
  return pyne::particle::is_heavy_ion(std::string(s));
}

bool pyne::particle::is_heavy_ion(std::string s) {
  if(pyne::nucname::isnuclide(s)) {
    if(pyne::particle::is_hydrogen(s))
      return false;
    else
      return true;
  }
  return false;
}

// is valid functions
bool pyne::particle::is_valid(int s) {
  if(pyne::nucname::isnuclide(s))
    return true;
  else
    return pyne::particle::is_valid(std::string(id_name[s]));
}

bool pyne::particle::is_valid(char *s) {
  return pyne::particle::is_valid(std::string(s));
}

bool pyne::particle::is_valid(std::string s) {
  // check std name
  if(0 < names.count(s))
    return true;
  // check alternative name
  if(0 < altnames.count(s))
    return true;
  // check if is a heavy ion
  if(pyne::nucname::isnuclide(s))
    return true;
  else
    return false;
}

// pdc functions
int pyne::particle::id(int s) {
  if (0 < pdc_nums.count(s))
    return s;
  else
    return 0;
}

int pyne::particle::id(char *s) {
  return pyne::particle::id(std::string(s));
}

int pyne::particle::id(std::string s) {
  if(pyne::nucname::isnuclide(s))
    {
      if(pyne::particle::is_hydrogen(s))
  return name_id["Proton"];
      if( pyne::particle::is_heavy_ion(s) )
  return 0;
    }

  if (0 < pdc_nums.count(name_id[s]))
    return name_id[s];
  if (0 < pdc_nums.count(altnames[s]))
    return altnames[s];
  return 0;
}

// name functions
std::string pyne::particle::name(int s) {
  if(s < 9999999)
    return pyne::particle::name(id_name[s]);
  if(pyne::nucname::isnuclide(s))
    return pyne::particle::name(pyne::nucname::name(s));
  return pyne::particle::name(id_name[s]);
}

std::string pyne::particle::name(char *s) {
  return pyne::particle::name(std::string(s));
}

std::string pyne::particle::name(std::string s) {
  // check if is a hydrogen
  if(pyne::nucname::isnuclide(s))
    {
      if(pyne::particle::is_hydrogen(s))
  return "Proton";
      if( pyne::particle::is_heavy_ion(s) )
  return s;
    }
  // check std name
  if(0 < names.count(s))
    return s;
  // check alternative name
  if(0 < altnames.count(s))
    return id_name[altnames[s]];
  // check for heavy ion
  else
    throw NotAParticle(s);
}

// convert name to mcnp id
std::string pyne::particle::mcnp(int s) {
  return pyne::particle::mcnp(pyne::particle::name(s));
}

std::string pyne::particle::mcnp(char *s) {
  return pyne::particle::mcnp(pyne::particle::name(s));
}

std::string pyne::particle::mcnp(std::string s) {
  if(0 < part_to_mcnp.count(pyne::particle::name(s)))
    return part_to_mcnp[pyne::particle::name(s)];
  else
    {
      std::cout << "Not a valid MCNP5 particle" << std::endl;
      return "?";
    }
}

// convert name to mcnp6 id
std::string pyne::particle::mcnp6(int s) {
  return pyne::particle::mcnp6(pyne::particle::name(s));
}

std::string pyne::particle::mcnp6(char *s) {
  return pyne::particle::mcnp6(pyne::particle::name(s));
}

std::string pyne::particle::mcnp6(std::string s) {
  if(0 < part_to_mcnp6.count(pyne::particle::name(s)))
    return part_to_mcnp6[pyne::particle::name(s)];
  else
    {
      std::cout << "Not a valid MCNP6 particle" << std::endl;
      return "?";
    }
}

// convert name to fluka id
std::string pyne::particle::fluka(int s) {
  return pyne::particle::fluka(pyne::particle::name(s));
}

std::string pyne::particle::fluka(char *s) {
  return pyne::particle::fluka(pyne::particle::name(s));
}

std::string pyne::particle::fluka(std::string s) {
  if (pyne::particle::is_heavy_ion(s))
    return "HEAVYION";
  else if(0 < part_to_fluka.count(pyne::particle::name(s)))
    return part_to_fluka[pyne::particle::name(s)];
  else
    {
      std::cout << "Not a valid Fluka particle" << std::endl;
      return "???????";
    }
}

// convert name to geant4 id
std::string pyne::particle::geant4(int s) {
  return pyne::particle::geant4(pyne::particle::name(s));
}

std::string pyne::particle::geant4(char *s) {
  return pyne::particle::geant4(pyne::particle::name(s));
}

std::string pyne::particle::geant4(std::string s) {
  if (pyne::particle::is_heavy_ion(s))
    return "GenericIon";
  else if(0 < part_to_geant4.count(pyne::particle::name(s)))
    return part_to_geant4[pyne::particle::name(s)];
  else
    {
      std::cout << "Not a valid Geant4 particle" << std::endl;
      return "???????";
    }
}


// describe functions
std::string pyne::particle::describe(int s) {
  if(pyne::nucname::isnuclide(s))
    return pyne::particle::describe(pyne::nucname::name(s));
  return pyne::particle::describe(id_name[s]);
}

std::string pyne::particle::describe(char *s) {
  return pyne::particle::describe(std::string(s));
}

std::string pyne::particle::describe(std::string s) {
  // check if is a hydrogen
  if (pyne::nucname::isnuclide(s))
    {
      if (pyne::particle::is_hydrogen(s))
  return docs[pyne::particle::name(s)];
      if (pyne::particle::is_heavy_ion(s))
  return "Is a heavy ion";
    }
  // check std name
  if(0 < names.count(s))
    return docs[s];
  // check alternative name
  if(0 < altnames.count(s))
    return docs[id_name[altnames[s]]];
  // check if is a heavy ion
  else
    throw NotAParticle(s);
}
//
// end of /home/mouginot/work/app/pyne/src/particle.cpp
//


//
// start of /home/mouginot/work/app/pyne/src/tally.cpp
//
// Tally.cpp
// Central Tally Class
// -- Andrew Davis

#include <iomanip>
#include <string>
#include <vector>

#ifndef PYNE_IS_AMALGAMATED
  #include "particle.h"
  #include "tally.h"
#endif

enum entity_type_enum {VOLUME, SURFACE, MESH}; // Enumeration for entity types
enum tally_type_enum {FLUX, CURRENT};  // Enumeration for tally types

const std::string tally_type_enum2string[] = {"Flux", "Current"};
const std::string entity_type_enum2string[] = {"Volume", "Surface", "Mesh"};
const std::string geometry_type_enum2string[] = {"Cartesian", "Cylinder"};

/***************************/
/*** Protected Functions ***/
/***************************/

// there are no protected functions currently
// fool.

/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/
pyne::Tally::Tally() {
  // Empty Tally Constructor
  tally_type = "";
  entity_id = -1;
  entity_type = "";
  entity_name = "";
  tally_name = "";
  entity_size = -1.0;
  normalization = 1.0;
}

// Default constructor
pyne::Tally::Tally(std::string type, std::string part_name, 
       int ent, std::string ent_type, 
       std::string ent_name, std::string tal_name,
       double size, double norm ) {

  tally_type = type;
  particle_names.push_back(pyne::particle::name(part_name));
  entity_id = ent;
  entity_type = ent_type;
  entity_name = ent_name;
  tally_name = tal_name;
  entity_size = size;
  normalization = norm;
}

pyne::Tally::Tally(std::string type, std::vector<std::string> part_names, 
       int ent, std::string ent_type, 
       std::string ent_name, std::string tal_name,
       double size, double norm ) {
  
  tally_type = type;
  particle_names = part_names;
  for (int i = 0; i < particle_names.size(); i++){
    particle_names[i] = pyne::particle::name(particle_names[i]);
  }
  entity_id = ent;
  entity_type = ent_type;
  entity_name = ent_name;
  tally_name = tal_name;
  entity_size = size;
  normalization = norm;
}

pyne::Tally::Tally(std::string part_name, std::string ent_geom,
                   std::vector<double> orgn,
                   std::vector<double> mesh_i, std::vector<double> mesh_j, std::vector<double> mesh_k,
                   std::vector<int> ints_i, std::vector<int> ints_j, std::vector<int> ints_k,
                   std::vector<double> e_bounds_, std::vector<int> e_ints_,
                   std::vector<double> axs_, std::vector<double> vec_,
                   std::string tal_name, double norm) {
  // Empty Tally Constructor
  entity_type = "Mesh";
  entity_name = "";
  particle_names.push_back(pyne::particle::name(part_name));
  entity_geometry = ent_geom;
  tally_name = tal_name;
  entity_size = -1;

  origin = orgn;
  vec = vec_;
  axs = axs_;
  meshes[0] = mesh_i;
  meshes[1] = mesh_j;
  meshes[2] = mesh_k;
  ints[0] = ints_i;
  ints[1] = ints_j;
  ints[2] = ints_k;
  e_bounds = e_bounds_;
  e_ints = e_ints_;
  normalization = norm;
}

// Destructor
pyne::Tally::~Tally() {}

/*--- Method definitions ---*/
//
void pyne::Tally::from_hdf5(char * filename, char *datapath, int row) {
  std::string fname(filename);
  std::string dpath(datapath);
  from_hdf5(fname,dpath,row);
}

//
void pyne::Tally::from_hdf5(std::string filename, std::string datapath, 
          int row) { 
  // line of data to acces
  int data_row = row;

  // check for file existence
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  // check to make sure is a HDF5 file
  bool is_h5 = H5Fis_hdf5(filename.c_str());
  if (!is_h5)
    throw h5wrap::FileNotHDF5(filename);

  // Open file and dataset.
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t dset = H5Dopen2(file, datapath.c_str(), H5P_DEFAULT);

  // Get dataspace and allocate memory for read buffer.
  hid_t space = H5Dget_space(dset);
  int rank  = H5Sget_simple_extent_ndims(space);
  hsize_t dims[1]; // for length of dataset 

  // get the length of the dataset
  int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
  
  // determine if chunked
  hid_t prop = H5Dget_create_plist(dset);
  
  hsize_t chunk_dimsr[1];
  int rank_chunk;

  if(H5D_CHUNKED == H5Pget_layout(prop))
    rank_chunk = H5Pget_chunk(prop, rank, chunk_dimsr);
  
  // allocate memory for data from file
  tally_struct* read_data = new tally_struct[dims[0]];

  // if row number is larger than data set only give last element
  if ( row >= dims[0] )
    data_row = dims[0]-1;
    
  // Create variable-length string datatype.
  hid_t strtype = H5Tcopy(H5T_C_S1);
  int status  = H5Tset_size(strtype, H5T_VARIABLE);
  
  // Create the compound datatype for memory.
  hid_t memtype = create_memtype();
  
  // Create the compound datatype for the file
  hid_t filetype = create_filetype();
  
  // Read the data.
  status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_data);

  // unpack the data and set values
  entity_id = read_data[data_row].entity_id;
  entity_type = entity_type_enum2string[read_data[data_row].entity_type];
  tally_type = tally_type_enum2string[read_data[data_row].tally_type];
  particle_names = pyne::split_string(read_data[data_row].particle_name, ",");
  tally_name = std::string(read_data[data_row].tally_name);
  entity_name = std::string(read_data[data_row].entity_name);
  entity_size = read_data[data_row].entity_size;
  normalization = read_data[data_row].normalization;

  // close the data sets
  status = H5Dclose(dset);
  status = H5Sclose(space);
  status = H5Tclose(filetype);
  status = H5Fclose(file);

  // tidy up
  delete[] read_data;
 
}

// Dummy Wrapper around C Style Functions
void pyne::Tally::write_hdf5(char * filename, char * datapath) {
  std::string fname(filename);
  std::string groupname(datapath);
  write_hdf5(fname,groupname);
}

// create filetype
hid_t pyne::Tally::create_filetype() {
  herr_t status;  // iostatus

  // create string type
  hid_t strtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(strtype, H5T_VARIABLE);

  hid_t filetype = H5Tcreate(H5T_COMPOUND, 8 + 8 + 8 + 
           (3*sizeof(hvl_t)) + 8 + 8);
  status = H5Tinsert(filetype, "entity_id", 0, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "entity_type", 8, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "tally_type", 8 + 8, H5T_STD_I64BE);
  status = H5Tinsert(filetype, "particle_name", 8 + 8 + 8, strtype);
  status = H5Tinsert(filetype, "entity_name", 8 + 8 + 8 + 
         sizeof(hvl_t), strtype);
  status = H5Tinsert(filetype, "tally_name", 8 + 8 + 8 + 
         (2*sizeof(hvl_t)) , strtype);
  status = H5Tinsert(filetype, "entity_size", 8 + 8 + 8 + 
         (3*sizeof(hvl_t)), H5T_IEEE_F64BE);
  status = H5Tinsert(filetype, "normalization", 8 + 8 + 8 + 
         (3*sizeof(hvl_t)) + 8, H5T_IEEE_F64BE);
  return filetype;
}

// create memory type 
hid_t pyne::Tally::create_memtype() {
  // iostatus
  herr_t status;

  //
  hid_t strtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(strtype, H5T_VARIABLE);
  // Create the compound datatype for memory.
  hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(tally_struct));
  status = H5Tinsert(memtype, "entity_id",
         HOFFSET(tally_struct, entity_id), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "entity_type",
         HOFFSET(tally_struct, entity_type), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "tally_type",
         HOFFSET(tally_struct, tally_type), H5T_NATIVE_INT);
  status = H5Tinsert(memtype, "particle_name",
         HOFFSET(tally_struct, particle_name), strtype);
  status = H5Tinsert(memtype, "entity_name",HOFFSET(tally_struct, entity_name),
         strtype);
  status = H5Tinsert(memtype, "tally_name",HOFFSET(tally_struct, tally_name),
         strtype);
  status = H5Tinsert(memtype, "entity_size",
         HOFFSET(tally_struct, entity_size), H5T_NATIVE_DOUBLE);
  status = H5Tinsert(memtype, "normalization",
         HOFFSET(tally_struct, normalization), H5T_NATIVE_DOUBLE);
  return memtype;
}

hid_t pyne::Tally::create_dataspace(hid_t file, std::string datapath) {
    // enable chunking 
    hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
    // set chunk size
    hsize_t chunk_dimensions[1]={1};
    herr_t status = H5Pset_chunk(prop, 1, chunk_dimensions);
    
    // allow varaible length strings
    hid_t strtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(strtype, H5T_VARIABLE);

    // Create the compound datatype for memory.
    hid_t memtype = create_memtype();
    
    // Create the compound datatype for the file
    hid_t filetype = create_filetype();
    
    // max dims unlimted
    hsize_t max_dims[1] = {H5S_UNLIMITED};
    // only ever let 1 tally object be added
    hsize_t dims[1] = {1}; 
    // Create dataspace.  Setting maximum size to NULL sets the maximum
    hid_t space = H5Screate_simple(1, dims, max_dims);

    // Create the dataset and write the compound data to it.
    return H5Dcreate2(file, datapath.c_str(), filetype, space, H5P_DEFAULT, prop,
           H5P_DEFAULT);
}

// Appends Tally object to dataset if file & datapath already exists
// if file exists & data path doesnt creates new datapath, 
// otherwise creates new file
void pyne::Tally::write_hdf5(std::string filename, std::string datapath) {

  // turn of annoying hdf5 errors
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  tally_struct tally_data[1]; // storage for the tally to add

  // setup the data to write
  tally_data[0].entity_id = entity_id;
  // entity type
  if (entity_type.find("Volume") != std::string::npos)
    tally_data[0].entity_type = VOLUME;
  else if (entity_type.find("Surface") != std::string::npos)
    tally_data[0].entity_type = SURFACE;

  // tally kind
  if (tally_type.find("Flux") != std::string::npos)
    tally_data[0].tally_type = FLUX;
  else if (tally_type.find("Current") != std::string::npos)
    tally_data[0].tally_type = CURRENT;

  // unpack from class to struct array
  std::string particle_name = pyne::join_to_string(particle_names, ",");
  tally_data[0].entity_id = entity_id;
  tally_data[0].entity_name = entity_name.c_str();
  tally_data[0].particle_name = particle_name.c_str();
  tally_data[0].tally_name = tally_name.c_str();
  tally_data[0].entity_size = entity_size;
  tally_data[0].normalization = normalization;

  
  // check for file existence
  bool is_exist = pyne::file_exists(filename);
  // create new file

  // check to make sure is a HDF5 file
  bool is_h5 = H5Fis_hdf5(filename.c_str());

  if (is_exist && !is_h5)
    throw h5wrap::FileNotHDF5(filename);

  if (!is_exist ) { // is a new file        
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, 
         H5P_DEFAULT);
    // create a dataspace
    hid_t dset = create_dataspace(file, datapath);

    hid_t memtype = create_memtype();
  
    herr_t status; // iostatus

    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tally_data);

    // close the data sets
    status = H5Dclose(dset);
    //    status = H5Sclose(space);
    //    status = H5Tclose(filetype);
    status = H5Fclose(file);
  
  }  else if ( is_exist && is_h5 ) {// already exists and is an hdf file
    // then we append the data to the end
    herr_t data_status; // iostatus

    // Open file and dataset.
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    // see if path exists
    data_status = H5Gget_objinfo (file, datapath.c_str(), 0, NULL);

    hid_t dset;
    // if fails neet to create dataset
    // still need to check that the datapath exists
    if (data_status != 0) // doesnt exist
      {
  dset = create_dataspace(file,datapath.c_str());
  hid_t memtype = create_memtype();
  herr_t status; // iostatus

  status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tally_data);

  // close the data sets
  status = H5Dclose(dset);
  status = H5Fclose(file);
      } else {
      
        dset = H5Dopen2(file, datapath.c_str(), H5P_DEFAULT);
      
  // Get dataspace and allocate memory for read buffer.
  hid_t space = H5Dget_space(dset);
  int rank  = H5Sget_simple_extent_ndims(space);
  hsize_t dims[1]; // for length of dataset 
  
  // get the length of the dataset
  int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
  
  // determine if chunked
  hid_t prop = H5Dget_create_plist(dset);

  hsize_t chunk_dimsr[1];
  int rank_chunk;
  if (H5D_CHUNKED == H5Pget_layout(prop))
    rank_chunk = H5Pget_chunk(prop, rank, chunk_dimsr);
  
  // allocate memory for data from file
  tally_struct* read_data = new tally_struct[dims[0]];
  
  // Create variable-length string datatype.
  hid_t strtype = H5Tcopy(H5T_C_S1);
  int status  = H5Tset_size(strtype, H5T_VARIABLE);
  
  // Create the compound datatype for memory.
  hid_t memtype = create_memtype();
  
  // Create the compound datatype for the file
  hid_t filetype = create_filetype();
  
  // Read the data.
  status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_data);
      
  // resize dims
  dims[0] += 1;
  
  // Extend the dataset
  status = H5Dextend(dset,dims);
  hid_t filespace = H5Dget_space(dset);
  // calculate the existing offset
  hsize_t offset[1] = {dims[0] - 1};  
  
  // select hyerslab
  hsize_t new_length[1] = {1};
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,offset , NULL,
             new_length, NULL);
  
  // create dataspace for new data
  space = H5Screate_simple(1,new_length, NULL);
  
  // Write the dataset to memory
  status = H5Dwrite(dset, memtype, space, filespace, H5P_DEFAULT, tally_data);
  
  // tidy up
  status = H5Dvlen_reclaim(memtype, space, H5P_DEFAULT, read_data);
  delete[] read_data;
  status = H5Dclose(dset);
  status = H5Sclose(space);
  status = H5Tclose(memtype);
  status = H5Tclose(strtype);
  status = H5Fclose(file);
    }
  }
}

std::ostream& operator<<(std::ostream& os, pyne::Tally tal) {
  //print the Tally to ostream
  os << "\t---------\n";
  os << "\t Tallying " << pyne::join_to_string(tal.particle_names, ", ") << " " << tal.tally_type << std::endl;
  os << "\t in/on " << tal.entity_type << " " << tal.entity_id << std::endl;
  return os;
}

// Sets string to valid mcnp formatted tally
// Takes mcnp version as arg, like 5 or 6
std::string pyne::Tally::mcnp(int tally_index, std::string mcnp_version,
                              std::string out) {
  std::stringstream output;  // output stream
  std::string particle_token = "";

  if (particle_names.size() == 0) {
      particle_token = "?";
  }

  // particle token
  for (int i = 0; i < particle_names.size(); i++ ){
    if (i > 0)
      particle_token += ",";
    if (mcnp_version.find("mcnp5") != std::string::npos)
      particle_token += pyne::particle::mcnp(particle_names[i]);
    else if (mcnp_version.find("mcnp6") != std::string::npos)
      particle_token += pyne::particle::mcnp6(particle_names[i]);
    else
      particle_token += "?";
  }
  // print out comment line
  output << "C " << tally_name << std::endl;
  output << std::setiosflags(std::ios::fixed) << std::setprecision(6);

  if (normalization != 1.0) output << std::scientific;
  int tally_id = 0;

  // neednt check entity type
  if (entity_type.find("Surface") != std::string::npos) {
    if (tally_type.find("Current") != std::string::npos) {
      tally_id = 1;
    } else if (tally_type.find("Flux") != std::string::npos) {
      tally_id = 2;
    }
    output << form_mcnp_tally(tally_index, tally_id, particle_token, entity_id,
                              entity_size, normalization);

  } else if (entity_type.find("Volume") != std::string::npos) {
    if (tally_type.find("Flux") != std::string::npos) {
      tally_id = 4;
    } else if (tally_type.find("Current") != std::string::npos) {
      // makes no sense in mcnp
      return "";
    }
    output << form_mcnp_tally(tally_index, tally_id, particle_token, entity_id,
                              entity_size, normalization);

  } else if (entity_type.find("Mesh") != std::string::npos) {
    output << form_mcnp_meshtally(tally_index, particle_token, entity_geometry,
                                  axs, vec, origin, meshes, ints, e_bounds,
                                  e_ints, out);

  } else {
    std::cout << "tally/entity combination makes no sense for MCNP"
              << std::endl;
  }
  // print sd card if area/volume specified
  return output.str();
}

template <typename T>
bool pyne::Tally::is_zero(T vect) {
  int size = sizeof(vect) / sizeof(vect[0]);
  bool result = true;
  for (int i = 0; i < size; i++) result &= (vect[i] == 0);
  return result;
}


// Form the tally line as function of its properties
std::string pyne::Tally::form_mcnp_tally(int tally_index, int type,
                                         std::string particle_token,
                                         int entity_id, double entity_size,
                                         double normalization) {
  std::stringstream tally_stream;  // tally stream
  tally_stream << std::setiosflags(std::ios::fixed) << std::setprecision(6);
  if (normalization != 1.0) tally_stream << std::scientific;

  tally_stream << "F" << tally_index << type << ":" << particle_token << " "
               << entity_id << std::endl;

  if (entity_size > 0.0)
    tally_stream << "SD" << tally_index << type << " " << entity_size
                 << std::endl;

  if (normalization != 1.0)
    tally_stream << "FM" << tally_index << type << " " << normalization
                 << std::endl;

  return tally_stream.str();
}


// Form the mesh tally line as function of its properties
std::string pyne::Tally::form_mcnp_meshtally(
    int tally_index, std::string particle_token, std::string entity_geometry,
    std::vector<double> axs, std::vector<double> vec,
    std::vector<double> origin, std::vector<double> meshes[3],
    std::vector<int> ints[3], std::vector<double> e_bounds,
    std::vector<int> e_ints, std::string out) {
  std::stringstream mtally_stream;
  // indentation block
  std::string indent_block = "          ";

  mtally_stream << "FMESH" << tally_index << "4:" << particle_token << " ";
  mtally_stream << "GEOM=";

  if (entity_geometry.find("Cartesian") != std::string::npos) {
    mtally_stream << "XYZ" << std::endl;
    mtally_stream << indent_block;
  } else if (entity_geometry.find("Cylinder") != std::string::npos) {
    mtally_stream << "CYL" << std::endl;
    if (!is_zero(axs)) {
      mtally_stream << indent_block << "AXS=" << pyne::join_to_string(axs) << std::endl;
    }
    if (!is_zero(vec)) {
      mtally_stream << indent_block << "VEC=" << pyne::join_to_string(vec) << std::endl;
    }
    mtally_stream << indent_block;
  }

  mtally_stream << "ORIGIN=" << pyne::join_to_string(origin) << std::endl;
  std::string dir_name[3] = {"I", "J", "K"};

  for (int j = 0; j < 3; j++) {
    mtally_stream << indent_block;
    mtally_stream << dir_name[j] << "MESH=" << pyne::join_to_string(meshes[j]);
    if (ints[j].size() > 0) {
      mtally_stream << " " << dir_name[j]
                    << "INTS=" << pyne::join_to_string(ints[j]);
    }
    mtally_stream << std::endl;
  }
  if (e_bounds.size() > 0) {
    mtally_stream << indent_block << "EMESH=" << pyne::join_to_string(e_bounds);
  }
  mtally_stream << std::endl;
  if (e_ints.size() > 0) {
    mtally_stream << indent_block << "EINTS=" << pyne::join_to_string(e_ints);
  }
  if (out.size() > 0) {
    mtally_stream << std::endl << indent_block << "OUT=" << out;
  }
  return mtally_stream.str();
}


// Produces valid fluka tally
std::string pyne::Tally::fluka(std::string unit_number) {
  std::stringstream output; // output stream
  
  if (particle_names.size() != 1) {
    std::cout << "Fluka multi-particles tally have not been yet implemented!" << std:: endl;
    exit(1);
  }
  std::string particle_name = particle_names[0];

  // check entity type
  if (entity_type.find("Volume") != std::string::npos) {
    // ok
  }  else if (entity_type.find("Surface") != std::string::npos) {
      std::cout << "Surface tally not valid in FLUKA" << std::endl;
  } else {
      std::cout << "Unknown entity type" << std::endl;
  }

  output << "* " << tally_name << std::endl;
  output << std::setiosflags(std::ios::fixed) << std::setprecision(1);
  // check tally type
  if (tally_type.find("Flux") != std::string::npos) {
      output << std::setw(10) << std::left  << "USRTRACK";
      output << std::setw(10) << std::right << "     1.0";
      output << std::setw(10) << std::right 
             << pyne::particle::fluka(particle_name);
      output << std::setw(10) << std::right << unit_number;
      output << std::setw(10) << std::right << entity_name;
      if(entity_size > 0.0) {
        output << std::scientific;
        output << std::setprecision(4);
        output << std::setw(10) << std::right << entity_size;
      }
      else
        output << std::setw(10) << std::right << 1.0;

      output << std::setw(10) << std::right << "   1000."; // number of eints
      tally_name.resize(8);
      output << std::setw(8) << std::left 
             << tally_name; // may need to make sure less than 10 chars
      output << std::endl;
      output << std::setw(10) << std::left  << "USRTRACK";
      output << std::setw(10) << std::right << "   10.E1";
      output << std::setw(10) << std::right << "   1.E-3";
      output << std::setw(10) << std::right << "        ";
      output << std::setw(10) << std::right << "        ";
      output << std::setw(10) << std::right << "        ";
      output << std::setw(10) << std::right << "        ";
      output << std::setw(8) << std::left << "       &";      
      // end of usrtrack
  } else if (tally_type.find("Current") != std::string::npos) {
      output << std::setw(10) << std::left  << "USRBDX  ";    
      output << std::setw(10) << std::right << "   110.0";
      output << std::setw(10) << std::right 
             << pyne::particle::fluka(particle_name);
      output << std::setw(10) << std::right << unit_number;
      output << std::setw(10) << std::right << entity_name; // upstream
      output << std::setw(10) << std::right << entity_name; // downstream
      if ( entity_size > 0.0 )
        output << std::setw(10) << std::right << entity_size; // area
      else
        output << std::setw(10) << std::right << 1.0;

      tally_name.resize(8);
      output << std::setw(8) << std::right 
             << tally_name; // may need to make sure less than 10 chars
      output << std::endl;
      output << std::setw(10) << std::left  << "USRBDX  ";    
      output << std::setw(10) << std::right << "  10.0E1";
      output << std::setw(10) << std::right << "     0.0";
      output << std::setw(10) << std::right << "  1000.0"; // number of ints
      output << std::setw(10) << std::right << "12.56637"; // 4pi
      output << std::setw(10) << std::right << "     0.0";
      output << std::setw(10) << std::right 
             << "   240.0"; // number of angular ints
      output << std::setw(8) << std::left << "       &";      
      // end of usrbdx
  } else {
    std::cout << "Unknown tally type" << std::endl;
  }
  return output.str();
}
//
// end of /home/mouginot/work/app/pyne/src/tally.cpp
//


