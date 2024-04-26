CC        = cc
OPTFLAG1  = -O3
OPTFLAG2  = -lm
TARGET    = ../run/solc

.SUFFIXES:
.SUFFIXES: .o .c 

.c.o:
	$(CC) -c $(CFLAGS) $(OPTFLAG1) $(OPTFLAG2) $< -o $@

OBJS = \
	input.o \
	pointer_init.o \
	boundary_cell.o \
	cell_metrics.o \
	poi_gen.o \
	solver_PCG.o \
	outucd.o \
	allocate.o \
	main.o

HEADERS = \
	struct.h \
	struct_ext.h \
	pcg.h \
	pcg_ext.h \
	input.h \
	pointer_init.h \
	boundary_cell.h \
	cell_metrics.h \
	poi_gen.h \
	solver_PCG.h \
	allocate.h \
	outucd.h

all: $(TARGET)

$(TARGET): $(OBJS) 
	$(CC) $(CFLAGS) $(OPTFLAG1) -o $@ $(OBJS) $(OPTFLAG2)

$(OBJS): $(HEADERS)

clean:
	rm -f *.o $(TARGET) *.log *~ *.lst
