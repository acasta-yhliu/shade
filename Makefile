CC = clang++
CXXFLAGS = -O3 -fopenmp -Iinclude -D_CRT_SECURE_NO_WARNINGS
SOURCES  = 
TARGET   = Shade.exe
MAIN     = src/ShadeMain.cc
CXXDFLAGS = -O0 -g -Iinclude -D_CRT_SECURE_NO_WARNINGS
TEST     = src/ShadeMain.cc
TESTTARGET = ShadeTest.exe

all:
	$(CC) $(CXXFLAGS) $(SOURCES) $(MAIN) -o $(TARGET)

test:
	$(CC) $(CXXDFLAGS) $(SOURCES) $(TEST) -o $(TESTTARGET)

run:
	make
	./Shade
	make view

view:
	irfanview image.ppm