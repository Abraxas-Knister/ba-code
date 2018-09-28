CC = g++
include .positionofeigen.mk
OP = -I $(Eigen) -I $(PWD)/src -Wall -Wextra -pedantic-errors -std=c++14 -O3 -march=native -fopenmp
SRC = $(shell ls src/*.cpp) 
OBJ = $(SRC:%.cpp=%.o)
MN = $(shell ls main/*.hpp)
DEPFILE = .depend

.PHONY: all clean dep
.PRECIOUS: %.o

all: $(MN:%.hpp=%.prog)

clean:
	-rm main/*.prog
	-rm src/*.o

dep:
	-$(CC) -MM -I $(PWD)/src $(SRC) > $(DEPFILE)
	#-$(CC) -MM -I $(PWD)/src $(MN) >> $(DEPFILE)

-include $(DEPFILE)
%.o: %.cpp
	$(CC) $(OP) -c $< 
	-mv *.o src/

%.prog: %.hpp $(OBJ)
	echo '#include "'$<'"' > main.cpp
	$(CC) $(OP) $(OBJ) main.cpp -o $@
	rm main.cpp
