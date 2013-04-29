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
    gzFile m_gzLogFileStream;
    bool m_bLogIsGz;
    string m_sLogFile;
    uint m_nIteration;

    // Wimpute redefinition of hmm_like
    // so far it only adds logging
    virtual  fast hmm_like(uint I, uint *P);

public:

    // print out usage
    static void document(void);
    
    // are we logging?
    bool LogOn() { return !m_sLogFile.empty(); }

    // set and open log file
    void SetLog( const string &sLogFile);
    
    // function for logging a string
    template <typename T>
    void WriteToLog (const T &tInput);

    void estimate(void);
};

// template definition of public function must go in header
template <typename T>
void Wimpute::WriteToLog( const T & tInput )
{

    stringstream ssTemp;
    ssTemp << tInput;
        
    if(m_bLogIsGz){
        gzprintf(m_gzLogFileStream, ssTemp.str().c_str());
    }
    else{
    // open logFileStream for append if it's not yet open
    if( !m_ofsLogFileStream ){
        m_ofsLogFileStream.open(m_sLogFile);
    }

    // exit if the file cannot be opened
    if( !m_ofsLogFileStream.is_open() ){
        cerr << "could not open log file "<< m_sLogFile <<" for writing" << endl;
        exit(1);
    }
    
    // write input to log file
    m_ofsLogFileStream << ssTemp.str();
    m_ofsLogFileStream.flush();
    }
};


#endif /* _WIMPUTE_H */

