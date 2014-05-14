# INSTALLATION
To compile all code use

    make

Always clean up directory before using make with a new option

    make clean

To compile without cuda support

    make NCUDA=1

To compile with debug lines and assertions turned on

    make debug

To compile and run all tests (also runs debug)

    make test
