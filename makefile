CC=g++

INCLUDES=-I/media/ssd/ngs-data-analysis/code/WFA/
ifneq ($(DEBUG),)
	CFLAGS=-g -fopenmp -Wall -pthread -std=c++14 -O0 -DDBGPRINT
else
	CFLAGS=-g -fopenmp -Wall -pthread -std=c++14 -O3 -mavx2
endif

CC_LDFLAGS=-lz -ltbb -std=c++14

TARGETS=index accalign-cpu
CPUSRC=reference.cpp accalign.cpp ssw.c ssw_cpp.cpp embedding.cpp ksw2_extz2_sse.c
IDXSRC=index.cpp embedding.cpp
HEADERS=$(wildcard *.h)

all:	${TARGETS}

index:	${IDXSRC} ${HEADERS}
	${CC} ${IDXSRC} -o $@ ${CFLAGS} ${CC_LDFLAGS} ${INCLUDES} -L/media/ssd/ngs-data-analysis/code/WFA/build -lwfa

accalign-cpu: ${CPUSRC} ${HEADERS}
	${CC} ${CPUSRC} -o $@ ${CFLAGS} ${CC_LDFLAGS} ${INCLUDES} -L/media/ssd/ngs-data-analysis/code/WFA/build -lwfa

clean:
	rm ${TARGETS} *.o
