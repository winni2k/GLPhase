/* @(#)wimpute.h
 */

#ifndef _WIMPUTE_H
#define _WIMPUTE_H 1
/*
#include <string>
#include <fstream>
#include <iostream>
*/
#include "impute.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class Wimpute : public Impute
{
    
private:
    ofstream m_ofsLogFileStream;
    string m_sLogFile;
public:

    // print out usage
    static void document(void);
    
    // are we logging?
    bool LogOn() { return !m_sLogFile.empty(); }

    // set and open log file
    void SetLog( const string &sLogFile);
    
    // function for logging a string
    void WriteToLog (const string &sInput);

};

#endif /* _WIMPUTE_H */

