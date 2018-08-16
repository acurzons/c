

TARGET	=sedv2


CC	=gcc
CCFLAGS	= -O2 -Wall -g
LDFLAGS	=
LDLIBS	=-lgsl -lgslcblas -lm


SRC_DIR	=src
OBJ_DIR	=obj
BIN_DIR	=bin


SOURCES	=$(wildcard $(SRC_DIR)/*.c)
OBJECTS	=$(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SOURCES))

EXE	=$(BIN_DIR)/$(TARGET)

all	: $(EXE)

$(EXE)	: $(OBJECTS)
	$(CC) $(CCFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(OBJ_DIR)/%.o	: $(SRC_DIR)/%.c
	$(CC) $(CCFLAGS) -c $^ -o $@

