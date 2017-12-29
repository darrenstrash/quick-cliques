
BUILD_DIR = build
SRC_DIR   = src
BIN_DIR   = bin

CFLAGS = -Winline -O2 -std=c++11 -g
#CFLAGS = -Winline -DDEBUG_MESSAGE -O0 -std=c++0x -g

SOURCES_TMP += CliqueTools.cpp
SOURCES_TMP += MemoryManager.cpp
SOURCES_TMP += Algorithm.cpp
SOURCES_TMP += TomitaAlgorithm.cpp
SOURCES_TMP += AdjacencyListAlgorithm.cpp
SOURCES_TMP += HybridAlgorithm.cpp
SOURCES_TMP += DegeneracyAlgorithm.cpp
SOURCES_TMP += DegeneracyTools.cpp
SOURCES_TMP += Tools.cpp

SOURCES=$(addprefix $(SOURCES_DIR)/, $(SOURCES_TMP))

OBJECTS_TMP=$(SOURCES_TMP:.cpp=.o)
OBJECTS=$(addprefix $(BUILD_DIR)/, $(OBJECTS_TMP))

DEPFILES_TMP:=$(SOURCES_TMP:.cpp=.d)
DEPFILES=$(addprefix $(BUILD_DIR)/, $(DEPFILES_TMP))

EXEC_NAMES = printnm compdegen qc

EXECS = $(addprefix $(BIN_DIR)/, $(EXEC_NAMES))

#DEFINE += -DDEBUG       #for debugging
#DEFINE += -DMEMORY_DEBUG #for memory debugging.
#DEFINE += -DPRINT_CLIQUES_ONE_BY_ONE   #print cliques, one per line

# print cliques in tree-like format:
#  - print each vertex that's evaluated
#  - print c if you reach a maximal clique to report
#  - print b if you backtrack
#DEFINE += -DPRINT_CLIQUES_TOMITA_STYLE # used by Eppstein and Strash (2011)

#some systems handle malloc and calloc with 0 bytes strangely.
DEFINE += -DALLOW_ALLOC_ZERO_BYTES      # used by Eppstein and Strash (2011) 

#set CXX to generic g++ if not set.
CXX ?= g++

VPATH = src

.PHONY : all

all: $(EXECS)

.PHONY : clean

clean: 
	rm -rf $(EXECS) $(BUILD_DIR) $(BIN_DIR)

$(BIN_DIR)/printnm: printnm.cpp ${OBJECTS} makefile | ${BIN_DIR}
	$(CXX) $(CFLAGS) ${DEFINE} ${OBJECTS} $(SRC_DIR)/printnm.cpp -o $@

$(BIN_DIR)/compdegen: compdegen.cpp ${OBJECTS} makefile | ${BIN_DIR}
	$(CXX) $(CFLAGS) ${DEFINE} ${OBJECTS} $(SRC_DIR)/compdegen.cpp -o $@

$(BIN_DIR)/qc: main.cpp ${OBJECTS} makefile | ${BIN_DIR}
	$(CXX) $(CFLAGS) ${DEFINE} ${OBJECTS} $(SRC_DIR)/main.cpp -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h $(BUILD_DIR)/%.d makefile | $(BUILD_DIR)
	$(CXX) $(CFLAGS) ${DEFINE} -c $< -o $@

$(BUILD_DIR)/%.d: $(SRC_DIR)/%.cpp makefile | $(BUILD_DIR)
	$(CXX) $(CFLAGS) -MM -MT '$(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$<)' $< -MF $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

