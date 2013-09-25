CPPFLAGS = -Wall -O3
CC = g++-4.1
DFLAGS = -DSHR_REP

BIN_DIR = bin
SICTIN_SRC = /seqdata/people/robin/SICTIN/SICTIN/src
TARGETS = $(BIN_DIR)/modeller $(BIN_DIR)/runner $(BIN_DIR)/change

GSLFLAGS_C = -I/usr/include/gsl -I$(SICTIN_SRC)
GSLFLAGS_L = -L/usr/lib/gsl -lgsl -lgslcblas -lm

all: clean $(TARGETS)

clean: 
	-rm -rf $(TARGETS)
	-rm -rf $(BIN_DIR)/*.o

bed_utils.o: $(SICTIN_SRC)/bed_utils.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) -c $(SICTIN_SRC)/bed_utils.cpp -o $(BIN_DIR)/bed_utils.o

utils.o: $(SICTIN_SRC)/utils.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) -c $(SICTIN_SRC)/utils.cpp -o $(BIN_DIR)/utils.o

np_utils.o: src/np_utils.cpp $(SICTIN_SRC)/utils.h
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) -c src/np_utils.cpp -o $(BIN_DIR)/np_utils.o

zeb_utils.o: src/zeb_utils.cpp $(SICTIN_SRC)/utils.h
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) -c src/zeb_utils.cpp -o $(BIN_DIR)/zeb_utils.o

np_math.o: src/np_math.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) -c src/np_math.cpp -o $(BIN_DIR)/np_math.o

np_predictor.o: src/np_predictor.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) -c src/np_predictor.cpp -o $(BIN_DIR)/np_predictor.o

$(BIN_DIR)/np: utils.o bed_utils.o np_utils.o np_math.o np_predictor.o 
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o $(BIN_DIR)/np_utils.o $(BIN_DIR)/np_math.o $(BIN_DIR)/np_predictor.o \
src/np.cpp -o $(BIN_DIR)/np	

$(BIN_DIR)/modeller: utils.o bed_utils.o zeb_utils.o np_math.o
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o $(BIN_DIR)/zeb_utils.o $(BIN_DIR)/np_math.o \
src/modeller.cpp -o $(BIN_DIR)/modeller	

$(BIN_DIR)/runner: utils.o bed_utils.o zeb_utils.o np_math.o
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o $(BIN_DIR)/zeb_utils.o $(BIN_DIR)/np_math.o \
src/runner.cpp -o $(BIN_DIR)/runner	

$(BIN_DIR)/change: utils.o bed_utils.o zeb_utils.o np_math.o
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o $(BIN_DIR)/zeb_utils.o $(BIN_DIR)/np_math.o \
src/change.cpp -o $(BIN_DIR)/change	

$(BIN_DIR)/bin2bed: utils.o bed_utils.o zeb_utils.o np_math.o
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o $(BIN_DIR)/zeb_utils.o $(BIN_DIR)/np_math.o \
src/bin2bed.cpp -o $(BIN_DIR)/bin2bed	

$(BIN_DIR)/stretches: utils.o bed_utils.o zeb_utils.o np_math.o
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o $(BIN_DIR)/zeb_utils.o $(BIN_DIR)/np_math.o \
src/stretches.cpp -o $(BIN_DIR)/stretches	

$(BIN_DIR)/genRandData: 
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o \
src/generateRandomData.cpp -o $(BIN_DIR)/genRandData	

$(BIN_DIR)/checkBin: 
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o \
src/checkBinaries.cpp -o $(BIN_DIR)/checkBin	

$(BIN_DIR)/checkStart: 
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o \
src/checkStarts.cpp -o $(BIN_DIR)/checkStart	

$(BIN_DIR)/detectPairs: 
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o \
src/detectPairs.cpp -o $(BIN_DIR)/detectPairs	

$(BIN_DIR)/checkOddsFile: 
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L) \
$(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o \
src/checkOddsFile.cpp -o $(BIN_DIR)/checkOddsFile	


