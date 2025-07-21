# Makefile for Spinning Objects in C

# Detect OS
ifeq ($(OS),Windows_NT)
    EXE = spinning_objects.exe
    CC = gcc
    CFLAGS =
else
    EXE = spinning_objects
    CC = gcc
    CFLAGS = -lm
endif

all: $(EXE)

spinning_objects.exe: spinning_objects.c
	$(CC) spinning_objects.c -o spinning_objects.exe

spinning_objects: spinning_objects.c
	$(CC) spinning_objects.c -o spinning_objects $(CFLAGS)

clean:
	rm -f spinning_objects spinning_objects.exe 