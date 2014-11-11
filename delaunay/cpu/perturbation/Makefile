# Makefile to compile code using CGAL
# Set CGAL_DIR to the directory where the cgal includes and libraries are

CGAL_DIR=/usr/lib64

INCLUDES += -I/usr/share/include
LIBS += -L/usr/lib64 -lCGAL -lboost_thread-mt -lgmp -lmpfr -lm
FILES += genmtrand.c
CFLAGS += -O3 -frounding-math
CPP = g++

%: %.cpp
	${CPP} ${CFLAGS} ${INCLUDES} ${LIBS} ${FILES} -o $@ $<

