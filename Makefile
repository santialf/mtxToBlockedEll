all: blocked_ell

blocked_ell: blocked_ell.cpp
	g++ $(INC) blocked_ell.cpp -o execute $(LIBS)
	gcc mmio.h 

clean:
	rm -f execute

.PHONY: clean all test
