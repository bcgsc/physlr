SOURCES = $(wildcard *.cc)
HEADERS = $(wildcard *.h)

OBJECTS = $(patsubst &.cc,%.o,$(SOURCES))

CC = c++
CFLAGS = -std=c++0x -Wall -m64
LDFLAGS =

OLD_MODE = $(shell cat .buildmode 2> /dev/null)
ifeq ($(DEBUG),1)
CFLAGS += -g -O0
ifneq ($(OLD_MODE),debug)
$(shell echo debug >.buildmode)
endif
else
CFLAGS += -O3 -DNDEBUG
ifneq ($(OLD_MODE),release)
$(shell echo release >.buildmode)
endif
endif

all: physlr-indexlr

%.o: %.cpp $(HEADERS) .buildmode
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -o $@ -c $<

physlr-indexlr: $(OBJECTS) .buildmode
	$(CC) $(OBJECTS) $(CFLAGS) $(LDFLAGS) $(EXTRA_CFLAGS) -o $@


clean:
	$(RM) physlr-indexlr *.o .buildmode

.PHONY: all clean
