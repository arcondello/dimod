ROOT := ..
SRC := $(ROOT)/dimod/include/
CATCH2 := $(ROOT)/testscpp/Catch2/single_include/

all: catch2 test_main tests

tests: test_main.out
	./test_main 

test_main: test_main.cpp
	$(CXX) -std=c++11 -Wall -Werror -c test_main.cpp -I $(CATCH2) 
	$(CXX) -std=c++11 -Wall -Werror test_main.o tests/*.cpp -o test_main -I $(SRC) -I $(CATCH2) 

catch2:
	git submodule init
	git submodule update
