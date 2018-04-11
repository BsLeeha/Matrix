# compiler
CC = g++

# compiler flags
# -g adds debug information to the executable file
# -Wall turns on most, but not all, compiler warnings 
CFLAGS = -g -Wall -fexceptions -std=c++11

# target entry: "default" or "all"
default : entry

# create executable file lee: these object files needed
entry: main.o
	$(CC) $(CFLAGS) -o lee main.o 

# main.o generate
main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp


# make clean
# exe: executable file
# .o: object file
# .~: backup file
clean: 
	del *.exe *.o *.~





