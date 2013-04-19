//
//  main.cpp
//  wimpute
//
//  Created by winni on 19/04/2013.
//  Copyright (c) 2013 warren kretzschmar. All rights reserved.
//

// my first program in C++

#include <iostream>
#include "add.h" // this brings in the declaration for add()
using namespace std;

int add (int x, int y);

int main(int argc, const char * argv[])
{
    cout << "Hello World! " << add(1,2) << endl;

    cout << "bool:\t\t" << sizeof(bool) << " bytes" << endl;
    cout << "char:\t\t" << sizeof(char) << " bytes" << endl;
    cout << "wchar_t:\t" << sizeof(wchar_t) << " bytes" << endl;
    cout << "short:\t\t" << sizeof(short) << " bytes" << endl;
    cout << "int:\t\t" << sizeof(int) << " bytes" << endl;
    cout << "long:\t\t" << sizeof(long) << " bytes" << endl;
    cout << "float:\t\t" << sizeof(float) << " bytes" << endl;
    cout << "double:\t\t" << sizeof(double) << " bytes" << endl;
    cout << "long double:\t" << sizeof(long double) << " bytes" << endl;

    return 0;
}
