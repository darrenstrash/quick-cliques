
OBJECT_DIR = obj
SRC_DIR    = src
BIN_DIR    = bin

CFLAGS = -Winline -O0 -std=c++11 -g

SOURCES_TMP += DegeneracyIndependentSets.cpp
SOURCES_TMP += MaximumCliqueAlgorithm.cpp
SOURCES_TMP += IndependentSets.cpp
SOURCES_TMP += CacheEfficientDegeneracyVertexSets.cpp
SOURCES_TMP += DegeneracyVertexSets.cpp
SOURCES_TMP += AdjacencyListVertexSets.cpp
SOURCES_TMP += BronKerboschAlgorithm.cpp
SOURCES_TMP += MemoryManager.cpp
SOURCES_TMP += Algorithm.cpp
SOURCES_TMP += TomitaAlgorithm.cpp
SOURCES_TMP += AdjacencyListAlgorithm.cpp
SOURCES_TMP += TimeDelayAdjacencyListAlgorithm.cpp
SOURCES_TMP += TimeDelayMaxDegreeAlgorithm.cpp
SOURCES_TMP += TimeDelayDegeneracyAlgorithm.cpp
SOURCES_TMP += HybridAlgorithm.cpp
SOURCES_TMP += DegeneracyAlgorithm.cpp
SOURCES_TMP += FasterDegeneracyAlgorithm.cpp
SOURCES_TMP += DegeneracyTools.cpp
SOURCES_TMP += Tools.cpp

SOURCES=$(addprefix $(SOURCES_DIR)/, $(SOURCES_TMP))

OBJECTS_TMP=$(SOURCES_TMP:.cpp=.o)
OBJECTS=$(addprefix $(OBJECT_DIR)/, $(OBJECTS_TMP))

EXEC_NAMES = printnm compdegen qc

EXECS = $(addprefix $(BIN_DIR)/, $(EXEC_NAMES))

#DEFINE += -DDEBUG       #for debugging
#DEFINE += -DMEMORY_DEBUG #for memory debugging.
#DEFINE += -DRETURN_CLIQUES_ONE_BY_ONE 
#DEFINE += -DPRINT_CLIQUES_ONE_BY_ONE 

#DEFINE += -DPRINT_CLIQUES_TOMITA_STYLE # used by Eppstein and Strash (2011)

#some systems handle malloc and calloc with 0 bytes strangely.
DEFINE += -DALLOW_ALLOC_ZERO_BYTES# used by Eppstein and Strash (2011) 

VPATH = src

.PHONY : all

all: $(EXECS)

.PHONY : clean

clean: 
	rm -rf $(EXECS) $(OBJECT_DIR) $(BIN_DIR)

$(BIN_DIR)/printnm: printnm.cpp ${OBJECTS} | ${BIN_DIR}
	g++ ${DEFINE} ${OBJECTS} $(SRC_DIR)/printnm.cpp -o $@

$(BIN_DIR)/compdegen: compdegen.cpp ${OBJECTS} | ${BIN_DIR}
	g++ $(CFLAGS) ${DEFINE} ${OBJECTS} $(SRC_DIR)/compdegen.cpp -o $@

$(BIN_DIR)/qc: main.cpp ${OBJECTS} | ${BIN_DIR}
	g++ $(CFLAGS) ${DEFINE} ${OBJECTS} $(SRC_DIR)/main.cpp -o $@

$(OBJECTS): $(OBJECT_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJECT_DIR)
	g++ $(CFLAGS) ${DEFINE} -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $(OBJECT_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
