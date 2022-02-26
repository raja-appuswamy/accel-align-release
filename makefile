CC=g++

ifneq ($(DEBUG),)
	CFLAGS=-g -Wall -pthread -O0 -DDBGPRINT -isystem./WFA -L$./WFA/build -std=c++14
else
	CFLAGS=-g -Wall -pthread -O3 -isystem./WFA -mavx2 -L./WFA/build -std=c++14
endif

ACCLDFLAGS=./WFA/build/libwfa.a -lz -ltbb 
TARGETS=accindex accalign
CPUSRC=reference.cpp accalign.cpp embedding.cpp ksw2_extz2_sse.c
IDXSRC=index.cpp embedding.cpp
HEADERS=$(wildcard *.h) 

.PHONY: WFA all
all: WFA ${TARGETS}

WFA:
	$(MAKE) -C WFA clean all

accindex: ${IDXSRC} ${HEADERS}
	${CC} -o $@ ${IDXSRC} ${ACCLDFLAGS} ${CFLAGS} 

accalign: ${CPUSRC} ${HEADERS}
	${CC} -o $@ ${CPUSRC} ${ACCLDFLAGS} ${CFLAGS}

clean:
	$(MAKE) -C WFA clean
	rm ${TARGETS}
