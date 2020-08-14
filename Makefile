#CFLAGS = -g
#FFLAGS = -g
#F90FLAGS = -g



BDIR = bin/
ODIR = objects/

# GMT 4.5.18 meca objects
GMT = $(GMTHOME)/
GMT_INC = -I$(GMT)/include/ -I$(GMT)/src/meca/
GMT_LIBS = -L$(GMT)/lib/ -L$(NETCDFDIR)/lib/ -lpsl -lgmt -lnetcdf 

MECA_OBJS = $(GMT)/src/meca/utilmeca.o

CAT_OBJS = $(ODIR)/handle_catalog.o  $(ODIR)/handle_catalog_gmt.o $(ODIR)/michael_leasq.o $(MECA_OBJS)

INCLUDES = $(GMT_INC)

PROGS = $(BDIR)/merge_catalog $(BDIR)/bin_catalog $(BDIR)/m02dcfp $(BDIR)/calc_gr $(BDIR)/m02mag

all: dirs progs

progs:
	make $(PROGS)

dirs:
	mkdir -p $(BDIR) $(ODIR)

clean:
	rm $(ODIR)/*

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



$(ODIR)/%.o: %.c  $(HDR_FLS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

