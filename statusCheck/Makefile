NAME=main
INCDIR := include
SRCDIR := src

FLAGS = `${ROOTSYS}/bin/root-config --cflags --libs`
DEFINES  := -I. -I$(INCDIR)

SOURCES  += $(SRCDIR)/camera.cpp
SOURCES  += $(SRCDIR)/ReadFile.cpp

$(NAME): $(NAME).C $(SOURCES)
	g++ -g $(DEFINES) -mcmodel=medium -o $@ $^ ${FLAGS}

clean:
	rm -f *.o
