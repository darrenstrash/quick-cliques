
OBJECT_DIR = obj
SRC_DIR    = src
BIN_DIR    = bin

OBJECTS  = $(OBJECT_DIR)/LinkedList.o
OBJECTS += $(OBJECT_DIR)/MemoryManager.o
OBJECTS += $(OBJECT_DIR)/tomita_algorithm.o
OBJECTS += $(OBJECT_DIR)/adjlist_algorithm.o
OBJECTS += $(OBJECT_DIR)/hybrid_algorithm.o
OBJECTS += $(OBJECT_DIR)/degeneracy_algorithm.o
OBJECTS += $(OBJECT_DIR)/degeneracy_helper.o
OBJECTS += $(OBJECT_DIR)/misc.o

EXEC_NAMES = printnm compdegen tomita adjlist hybrid degeneracy

EXECS = $(addprefix $(BIN_DIR)/, $(EXEC_NAMES))

#DEFINE += -DDEBUG       #for debugging
#DEFINE += -DMEMORY_DEBUG #for memory debugging.
#DEFINE += -DRETURN_CLIQUES_ONE_BY_ONE 
#DEFINE += -DPRINT_CLIQUES_ONE_BY_ONE 

DEFINE += -DPRINT_CLIQUES_TOMITA_STYLE # used by Eppstein and Strash (2011)

#some systems handle malloc and calloc with 0 bytes strangely.
DEFINE += -DALLOW_ALLOC_ZERO_BYTES # used by Eppstein and Strash (2011) 

VPATH = src

.PHONY : all

all: $(EXECS)

.PHONY : clean

clean: 
	rm -rf $(OBJECTS) $(EXECS) $(OBJECT_DIR) $(BIN_DIR)

$(BIN_DIR)/printnm: printnm.c ${OBJECTS} ${BIN_DIR}
	gcc -O2 -g ${DEFINE} ${OBJECTS} $(SRC_DIR)/printnm.c -o $@

$(BIN_DIR)/compdegen: compdegen.c ${OBJECTS} ${BIN_DIR}
	gcc -O2 -g ${DEFINE} ${OBJECTS} $(SRC_DIR)/compdegen.c -o $@

$(BIN_DIR)/tomita: tomita.c ${OBJECTS} ${BIN_DIR}
	gcc -O2 -g ${DEFINE} ${OBJECTS} $(SRC_DIR)/tomita.c -o $@

$(BIN_DIR)/adjlist: adjlist.c ${OBJECTS} ${BIN_DIR}
	gcc -O2 -g ${DEFINE} ${OBJECTS} $(SRC_DIR)/adjlist.c -o $@

$(BIN_DIR)/hybrid: hybrid.c ${OBJECTS} ${BIN_DIR}
	gcc -O2 -g ${DEFINE} ${OBJECTS} $(SRC_DIR)/hybrid.c -o $@

$(BIN_DIR)/degeneracy: degeneracy.c ${OBJECTS} ${BIN_DIR}
	gcc -O2 -g ${DEFINE} ${OBJECTS} $(SRC_DIR)/degeneracy.c -o $@

$(OBJECT_DIR)/LinkedList.o: LinkedList.c LinkedList.h ${OBJECT_DIR}
	gcc -O2 -g ${DEFINE} -c $(SRC_DIR)/LinkedList.c -o $@

$(OBJECT_DIR)/MemoryManager.o: MemoryManager.c MemoryManager.h ${OBJECT_DIR}
	gcc -O2 -g ${DEFINE} -c $(SRC_DIR)/MemoryManager.c -o $@

$(OBJECT_DIR)/tomita_algorithm.o: tomita_algorithm.c tomita_algorithm.h ${OBJECT_DIR}
	gcc -O2 -g ${DEFINE} -c $(SRC_DIR)/tomita_algorithm.c -o $@

$(OBJECT_DIR)/adjlist_algorithm.o: adjlist_algorithm.c adjlist_algorithm.h ${OBJECT_DIR}
	gcc -O2 -g ${DEFINE} -c $(SRC_DIR)/adjlist_algorithm.c -o $@

$(OBJECT_DIR)/hybrid_algorithm.o: hybrid_algorithm.c hybrid_algorithm.h ${OBJECT_DIR}
	gcc -O2 -g ${DEFINE} -c $(SRC_DIR)/hybrid_algorithm.c -o $@

$(OBJECT_DIR)/degeneracy_algorithm.o: degeneracy_algorithm.c degeneracy_algorithm.h ${OBJECT_DIR}
	gcc -O2 -g ${DEFINE} -c $(SRC_DIR)/degeneracy_algorithm.c -o $@

$(OBJECT_DIR)/degeneracy_helper.o: degeneracy_helper.c degeneracy_helper.h ${OBJECT_DIR}
	gcc -O2 -g ${DEFINE} -c $(SRC_DIR)/degeneracy_helper.c -o $@

$(OBJECT_DIR)/misc.o: misc.c misc.h ${OBJECT_DIR}
	gcc -O2 -g ${DEFINE} -c $(SRC_DIR)/misc.c -o $@

 ${OBJECT_DIR}:
	mkdir ${OBJECT_DIR}

${BIN_DIR}:
	mkdir ${BIN_DIR}
