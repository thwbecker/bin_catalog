#CFLAGS = -g
#FFLAGS = -g
#F90FLAGS = -g
FFLAGS = -O2
CFLAGS = -O2

LDFLAGS = -lm
BDIR = bin/
ODIR = objects/

ARCH=$(shell uname -m | awk '{print(tolower($$1))}')

# for eigenroutines, my version of NETLIB eispack
EISLIB = -Leispack/$(ARCH)/ -lmyeis

# GMT 4.5.18 meca objects - used to be linked, now included
#GMT = $(GMT4HOME)/
#GMT_INC = -I$(GMT)/include/ -I$(GMT)/src/meca/
#GMT_LIBS = -L$(GMT)/lib/ -L$(NETCDFDIR)/lib/ -lpsl -lgmt -lnetcdf 

#MECA_OBJS =  $(GMT)/src/meca/utilmeca.o

CAT_OBJS = $(ODIR)/handle_catalog.o  $(ODIR)/handle_catalog_gmt.o $(ODIR)/michael_leasq.o $(MECA_OBJS)

#INCLUDES = $(GMT_INC) -DHAVE_SINCOS
INCLUDES = $(GMT_INC) 
#
# main programs
PROGS = $(BDIR)/merge_catalog $(BDIR)/bin_catalog $(BDIR)/m02dcfp $(BDIR)/calc_gr $(BDIR)/m02mag

# just needed for Simpson style stress state plotting
EIGEN_PROGS = $(BDIR)/eigen  $(BDIR)/eigenvalues $(BDIR)/eigen3ds  $(BDIR)/eigenvalues3ds 

all: dirs progs libs eigen_progs

progs:
	make $(PROGS)

eigen_progs:
	make $(EIGEN_PROGS)

libs:
	cd eispack; make ; cd ..

dirs:
	mkdir -p $(BDIR) $(ODIR)

clean:
	rm $(ODIR)/*

dist_clean:
	rm $(BDIR)/*

$(BDIR)/merge_catalog: merge_catalog.c $(CAT_OBJS)
	$(CC) $(CFLAGS) merge_catalog.c $(INCLUDES)  $(CAT_OBJS) \
	-o $(BDIR)/merge_catalog  $(GMT_LIBS) $(LDFLAGS)

$(BDIR)/bin_catalog: bin_catalog.c $(CAT_OBJS)  catalog.h
	$(CC) $(CFLAGS) bin_catalog.c $(INCLUDES)  $(CAT_OBJS) \
	-o $(BDIR)/bin_catalog    $(GMT_LIBS)  $(LDFLAGS)

$(BDIR)/m02mag: m02mag.c $(CAT_OBJS) catalog.h
	$(CC) $(CFLAGS) m02mag.c $(INCLUDES)   $(CAT_OBJS) \
	-o $(BDIR)/m02mag  $(GMT_LIBS)  $(LDFLAGS)

$(BDIR)/m02dcfp: m02dcfp.c $(CAT_OBJS)   catalog.h
	$(CC) $(CFLAGS) m02dcfp.c $(INCLUDES)  $(CAT_OBJS) \
	-o $(BDIR)/m02dcfp    $(GMT_LIBS)  $(LDFLAGS)

$(BDIR)/calc_gr: calc_gr.c  $(CAT_OBJS)   catalog.h
	$(CC) $(CFLAGS) calc_gr.c $(INCLUDES)  $(CAT_OBJS)  \
	-o $(BDIR)/calc_gr    $(GMT_LIBS)  $(LDFLAGS)

$(BDIR)/eigen: $(ODIR)/eigen.main.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/eigen.main.o $(ODIR)/eigen.o \
	-o $(BDIR)/eigen $(EISLIB) $(LDFLAGS) 

$(BDIR)/eigenvalues: $(ODIR)/eigen.ov.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/eigen.ov.o $(ODIR)/eigen.o -o $(BDIR)/eigenvalues \
	$(EISLIB) $(LDFLAGS) 

$(BDIR)/eigen3ds: $(ODIR)/eigen.tds.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/eigen.tds.o  $(ODIR)/eigen.o \
	-o $(BDIR)/eigen3ds $(EISLIB) $(LDFLAGS) 

$(BDIR)/eigenvalues3ds: $(ODIR)/eigen.tds.ov.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/eigen.tds.ov.o $(ODIR)/eigen.o \
	-o $(BDIR)/eigenvalues3ds \
	$(EISLIB) $(LDFLAGS) 

$(ODIR)/eigen.main.o: eigen_driver.c
	$(CC) $(CFLAGS) -c eigen_driver.c -o $(ODIR)/eigen.main.o 

$(ODIR)/eigen.ov.o: eigen_driver.c
	$(CC) $(CFLAGS) -c eigen_driver.c -DONLY_VALUES -o $(ODIR)/eigen.ov.o

$(ODIR)/eigen.tds.o: eigen_driver.c
	$(CC) $(CFLAGS) -c eigen_driver.c -DTHREED_SYMMETRIC -o $(ODIR)/eigen.tds.o

$(ODIR)/eigen.tds.ov.o: eigen_driver.c
	$(CC) $(CFLAGS) -c eigen_driver.c -DTHREED_SYMMETRIC \
		-DONLY_VALUES -o $(ODIR)/eigen.tds.ov.o


$(ODIR)/%.o: %.c  $(HDR_FLS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

