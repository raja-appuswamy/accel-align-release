CC=g++
WFA_PATH=/media/ssd/ngs-data-analysis/code/WFA/

INCLUDES=-I${WFA_PATH}
ifneq ($(DEBUG),)
	CFLAGS=-g -fopenmp -Wall -pthread -std=c++14 -O0 -DDBGPRINT
else
	CFLAGS=-g -fopenmp -Wall -pthread -std=c++14 -O3 -mavx2
endif

CC_LDFLAGS=-lz -ltbb -std=c++14

TARGETS=accindex accalign
CPUSRC=reference.cpp accalign.cpp embedding.cpp ksw2_extz2_sse.c
IDXSRC=index.cpp embedding.cpp
HEADERS=$(wildcard *.h)

all: ${TARGETS}

accindex: ${IDXSRC} ${HEADERS}
	${CC} ${IDXSRC} -o $@ ${CFLAGS} ${CC_LDFLAGS} ${INCLUDES} -L${WFA_PATH}/build -lwfa

accalign: ${CPUSRC} ${HEADERS}
	${CC} ${CPUSRC} -o $@ ${CFLAGS} ${CC_LDFLAGS} ${INCLUDES} -L${WFA_PATH}/build -lwfa

clean:
	rm ${TARGETS} *.o
