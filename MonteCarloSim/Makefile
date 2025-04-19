CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11

INCLUDES = -I. -I../external/freeglut/include
LIBS = -lGL -lGLU -lglut -lOpenCL

SOURCES = MonteCarloSim.cpp camera.cpp matrix4x4.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = MonteCarloSimApp

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o $(EXECUTABLE)
