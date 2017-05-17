include Makefile.inc

#####################################
## DIRECTORIES
#####################################
CLUSTER=cluster
DISTR=distr
EXEC=$(DISTR)/benchsolv
BUILD=build
SRC=src
OPTIONS=options
OBJECTS := $(patsubst %.cpp,%.o,$(wildcard $(PROJECT)/*.cpp))

#####################################
## RULES
#####################################
all: $(EXEC)

$(EXEC): $(BUILD)/Metrics.o $(BUILD)/Benchmark.o $(BUILD)/Solver.o $(BUILD)/Mumps.o $(BUILD)/QR_Mumps.o $(BUILD)/main.o
	mkdir -p $(DISTR)
	$(CC) -o $@ $^ $(LDFLAGS)

$(BUILD)/%.o: $(PROJECT)/$(SRC)/%.cpp
	mkdir -p $(BUILD)
	$(CC) $(CFLAGS) -MF "$@.d" -o $@ $<

clean:
	rm -rf $(BUILD)

cleaner: clean
	rm -rf $(EXEC)
