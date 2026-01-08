CXX = g++

CXXFLAGS = -O3 -march=native -ffast-math $(shell root-config --cflags)
LDFLAGS = $(shell root-config --libs)

HEADERS = LocalEquilibriumModel.h  SimulationTools.h

TARGET_DP = data_production.exe
SOURCE_DP = script_create_data.C

TARGET_PM = call_and_save.exe
SOURCE_PM = script_call_and_save.C

TARGET_DBL = do_deblurring.exe
SOURCE_DBL = script_deblurring.C

TARGET_ERR = error_calculation.exe
SOURCE_ERR = script_estimate_error.C

all: $(TARGET_DP) $(TARGET_PM) $(TARGET_DBL) $(TARGET_ERR)

$(TARGET_DP): $(SOURCE_DP) $(HEADERS)
	$(CXX) -o $@ $(SOURCE_DP) $(CXXFLAGS) $(LDFLAGS)

$(TARGET_PM): $(SOURCE_PM) $(HEADERS)
	$(CXX) -o $@ $(SOURCE_PM) $(CXXFLAGS) $(LDFLAGS)

$(TARGET_DBL): $(SOURCE_DBL) $(HEADERS)
	$(CXX) -o $@ $(SOURCE_DBL) $(CXXFLAGS) $(LDFLAGS)

$(TARGET_ERR): $(SOURCE_ERR) $(HEADERS)
	$(CXX) -o $@ $(SOURCE_ERR) $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f $(TARGET_DP) $(TARGET_PM) $(TARGET_DBL) $(TARGET_ERR)

.PHONY: all clean

