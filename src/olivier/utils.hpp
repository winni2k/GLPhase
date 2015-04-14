//$Id: utils.h 611 2012-07-18 14:14:00Z koskos $

#ifndef _UTILS_H
#define _UTILS_H

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define PI 3.14159265358979323846

#define LOW_POS_DOUBLE 1e-300
#define BIG_POS_DOUBLE 1e300
#define LOW_NEG_DOUBLE -1e-300
#define BIG_NEG_DOUBLE -1e300
#define LOW_POS_FLOAT 1e-8
#define BIG_POS_FLOAT 1e8
#define LOW_NEG_FLOAT -1e-8
#define BIG_NEG_FLOAT -1e8
#define BIG_POS_INT 1000000000
#define BIG_NEG_INT -1000000000

#include <string>
#include <vector>
#include <queue>
#include <map>
#include <bitset>
#include <list>
#include <unordered_map>
#include <bitset>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <exception>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>

// using namespace std;
namespace bio = boost::iostreams;
namespace bpo = boost::program_options;
namespace bid = boost::uuids;

/******************************************************/
/*                  UTILS STATISTICS                  */
/******************************************************/
namespace putils {
void initRandom(long s);
double getRandom();
std::string getRandomID();
int getRandom(int);
long getSeed();
void normalise(std::vector<double> &v);
int sample(std::vector<double> &v, double sum);
double entropy(std::vector<double> &v);
double KLdistance(std::vector<double> &P, std::vector<double> &Q);
};

/******************************************************/
/*                  UTILS ALGORITHM                   */
/******************************************************/
namespace autils {
int max(std::vector<double> &v);
int max(std::vector<int> &v);
void findUniqueSet(std::vector<bool> &B, std::vector<int> &U); // ?
void decompose(int min, std::vector<std::vector<int>> &B,
               std::vector<std::vector<std::vector<int>>> &BB); //?
int checkDuo(int pa1, int pa2, int ca1, int ca2);
int checkTrio(int fa1, int fa2, int ma1, int ma2, int ca1, int ca2);
};

/******************************************************/
/*                  UTILS STRING                      */
/******************************************************/
namespace sutils {
int tokenize(std::string &, std::vector<std::string> &);
int tokenize(std::string &, std::vector<std::string> &, int);
std::string uint2str(unsigned int n);
std::string int2str(int n);
std::string int2str(std::vector<int> &v);
std::string long2str(long int n);
std::string double2str(double n, int prc = 4);
std::string double2str(std::vector<double> &v, int prc = 4);
std::string bool2str(std::vector<bool> &v);
std::string date2str(time_t *t, std::string format);
};

/******************************************************/
/*                  UTILS FILE                        */
/******************************************************/
namespace futils {
bool isFile(std::string f);
bool createFile(std::string f);
std::string extensionFile(std::string &filename);
void bool2binary(std::vector<bool> &V, std::ostream &fd);
bool binary2bool(std::vector<bool> &V, std::istream &fd);
};

/******************************************************/
/*                  EXCEPTIONS                        */
/******************************************************/
class myException : public std::exception {
public:
  explicit myException(std::string msg) : msg_(msg) {}

  virtual ~myException() throw() {}

  virtual const char *what() const throw() { return msg_.c_str(); }

private:
  std::string msg_;
};

/******************************************************/
/*                  INPUT FILE                        */
/******************************************************/
class ifile : public bio::filtering_istream {
private:
  std::string m_ext;
  std::string file;
  std::ifstream fd;
  bool m_isGood = false;

public:
  ifile();
  ifile(std::string filename, bool binary = false, std::string ext = "");
  ~ifile();
  std::string name();
  bool open(std::string filename, bool binary = false, std::string ext = "");
  bool readString(std::string &);
  void close();
  bool isGood() { return m_isGood; };
};

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
class ofile : public bio::filtering_ostream {
private:
  std::string file;
  std::ofstream fd;

public:
  ofile();
  ofile(std::string filename, bool binary = false);
  ~ofile();
  std::string name();
  bool open(std::string filename, bool binary = false);
  void writeString(std::string &);
  void close();
};

/******************************************************/
/*                  LOG FILE                          */
/******************************************************/
class lfile {
private:
  std::string file;
  std::ofstream fd;
  bool verboseC;
  bool verboseL;

public:
  lfile();
  ~lfile();
  std::string name();
  bool open(std::string filename = "file.log");
  void close();
  std::string getPrefix();
  void muteL();
  void unmuteL();
  void muteC();
  void unmuteC();
  void print(std::string s);
  void printC(std::string s);
  void printL(std::string s);
  void println(std::string s);
  void printlnC(std::string s);
  void printlnL(std::string s);
  void warning(std::string s);
  void error(std::string s);
};

#endif
