
BUILD_DIR = build
SRC_DIR   = src
BIN_DIR   = bin

CFLAGS = -Winline -O2 -std=c++11 -g
#CFLAGS = -Winline -DDEBUG_MESSAGE -O0 -std=c++11 -g

#SOURCES_TMP += LightWeightSparseFullMCS.cpp
#SOURCES_TMP += LightWeightSparseStaticOrderMCS.cpp
#SOURCES_TMP += LightWeightSparseMCR.cpp

SOURCES_TMP += TesterMISS.cpp
SOURCES_TMP += TesterStaticOrderMISS.cpp
SOURCES_TMP += TesterMISQ.cpp
SOURCES_TMP += ConnectedComponentMISS2.cpp
SOURCES_TMP += ConnectedComponentMISS.cpp
SOURCES_TMP += LightWeightReductionDominationMISR.cpp
SOURCES_TMP += Reducer.cpp
SOURCES_TMP += LightWeightSparseMCQ.cpp
SOURCES_TMP += LightWeightReductionDominationMISQ.cpp
SOURCES_TMP += MaxSubgraphAlgorithm.cpp
SOURCES_TMP += IsolatesIndependentSetColoringStrategy.cpp
SOURCES_TMP += SparseIndependentSetColoringStrategy.cpp
SOURCES_TMP += LightWeightReductionSparseFullMISS.cpp
SOURCES_TMP += LightWeightReductionSparseStaticOrderMISS.cpp
SOURCES_TMP += LightWeightReductionSparseMISR.cpp
SOURCES_TMP += LightWeightReductionSparseMISQ.cpp
SOURCES_TMP += LightWeightReductionFullMISS.cpp
SOURCES_TMP += LightWeightReductionStaticOrderMISS.cpp
SOURCES_TMP += LightWeightReductionMISR.cpp
SOURCES_TMP += LightWeightReductionMISQ.cpp
SOURCES_TMP += LightWeightFullMISS.cpp
SOURCES_TMP += LightWeightStaticOrderMISS.cpp
SOURCES_TMP += LightWeightMISR.cpp
SOURCES_TMP += LightWeightMISQ.cpp
SOURCES_TMP += OrderingTools.cpp
SOURCES_TMP += LightWeightFullMCS.cpp
SOURCES_TMP += LightWeightStaticOrderMCS.cpp
SOURCES_TMP += LightWeightMCR.cpp
SOURCES_TMP += LightWeightMCQ.cpp
SOURCES_TMP += AdjacencyMatrixVertexSetsMax.cpp
SOURCES_TMP += AdjacencyMatrixVertexSets.cpp
SOURCES_TMP += AdjacencyListVertexSetsMax.cpp
SOURCES_TMP += SparseCliqueColoringStrategy.cpp
SOURCES_TMP += IndependentSetColoringStrategy.cpp
SOURCES_TMP += CliqueColoringStrategy.cpp
SOURCES_TMP += Isolates.cpp
SOURCES_TMP += Isolates2.cpp
SOURCES_TMP += Isolates3.cpp
SOURCES_TMP += Isolates4.cpp
SOURCES_TMP += IsolatesWithMatrix.cpp
SOURCES_TMP += ExperimentalReduction.cpp
SOURCES_TMP += IndependentSetsReduction.cpp
SOURCES_TMP += Staging.cpp
SOURCES_TMP += CliqueTools.cpp
SOURCES_TMP += GraphTools.cpp
SOURCES_TMP += DegeneracyIndependentSets2.cpp
SOURCES_TMP += MaximalCliqueAlgorithm.cpp
SOURCES_TMP += MinimumCliqueAlgorithm.cpp
SOURCES_TMP += ReverseDegeneracyVertexSets.cpp
SOURCES_TMP += CliqueGraphAlgorithm.cpp
SOURCES_TMP += PartialMatchGraph.cpp
SOURCES_TMP += PartialMatchDegeneracyVertexSets.cpp
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
#SOURCES_TMP += ComparisonFullMISS.cpp
#SOURCES_TMP += ComparisonStaticOrderMISS.cpp
#SOURCES_TMP += ComparisonMISQ.cpp

SOURCES=$(addprefix $(SOURCES_DIR)/, $(SOURCES_TMP))

OBJECTS_TMP=$(SOURCES_TMP:.cpp=.o)
OBJECTS=$(addprefix $(BUILD_DIR)/, $(OBJECTS_TMP))

DEPFILES_TMP:=$(SOURCES_TMP:.cpp=.d)
DEPFILES=$(addprefix $(BUILD_DIR)/, $(DEPFILES_TMP))

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
	rm -rf $(EXECS) $(BUILD_DIR) $(BIN_DIR)

$(BIN_DIR)/printnm: printnm.cpp ${OBJECTS} | ${BIN_DIR}
	g++ ${DEFINE} ${OBJECTS} $(SRC_DIR)/printnm.cpp -o $@

$(BIN_DIR)/compdegen: compdegen.cpp ${OBJECTS} | ${BIN_DIR}
	g++ $(CFLAGS) ${DEFINE} ${OBJECTS} $(SRC_DIR)/compdegen.cpp -o $@

$(BIN_DIR)/qc: main.cpp ${OBJECTS} | ${BIN_DIR}
	g++ $(CFLAGS) ${DEFINE} ${OBJECTS} $(SRC_DIR)/main.cpp -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h $(BUILD_DIR)/%.d | $(BUILD_DIR)
	g++ $(CFLAGS) ${DEFINE} -c $< -o $@

$(BUILD_DIR)/%.d: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	g++ $(CFLAGS) -MM -MT '$(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$<)' $< -MF $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

