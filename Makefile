#CFLAGS = -g
#DEBUG_DEF = -DDEBUG

CFLAGS = -O2
DEBUG_DEF =

F90FLAGS = $(CFLAGS)
FFLAGS =  $(CFLAGS)

# copy out if you do not have sincos
SINCOS_DEF = -DHAVE_SINCOS


DEFINES = $(SINCOS_DEF) $(DEBUG_DEF)

# cproto -f2 *.c | grep -v main > proto.h
HFILES = catalog.h proto.h eigen.h

LDFLAGS = -lm
BDIR = bin/
ODIR = objects/

ARCH=$(shell uname -m | awk '{print(tolower($$1))}')

# for eigenroutines, my version of NETLIB eispack
EISLIB = -Leispack/$(ARCH)/ -lmyeis



CAT_OBJS = $(ODIR)/handle_catalog.o  $(ODIR)/fault_eq.o   \
	$(ODIR)/linalg_misc_geo.o  $(ODIR)/bvalue.o \
	$(ODIR)/handle_catalog_gmt.o $(ODIR)/michael_leasq.o $(ODIR)/geo_kdtree.o \
	$(MECA_OBJS) 

SINV_OBS = 	$(ODIR)/stress_inversion.o $(ODIR)/stability_criterion.o $(ODIR)/slip_deviation.o  


INCLUDES = $(GMT_INC)
#INCLUDES = $(GMT_INC) 
#
# main programs
PROGS = $(BDIR)/merge_catalog $(BDIR)/bin_catalog $(BDIR)/nsample_catalog $(BDIR)/solve_stress_one_bin \
	$(BDIR)/m02dcfp $(BDIR)/calc_gr $(BDIR)/calc_gr_time $(BDIR)/m02mag

# just needed for Simpson style stress state plotting
EIGEN_PROGS = $(BDIR)/eigen  $(BDIR)/eigenvalues $(BDIR)/eigen3ds  $(BDIR)/eigenvalues3ds 

#
TEST_PROGS = $(BDIR)/test_eigen $(BDIR)/test_kdtree

all: dirs progs libs eigen_progs test_progs

progs:
	make $(PROGS)

eigen_progs:
	make $(EIGEN_PROGS)

test_progs:
	make $(TEST_PROGS)

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
	-o $(BDIR)/merge_catalog   $(LDFLAGS)


$(BDIR)/bin_catalog: bin_catalog.c $(CAT_OBJS)  $(SINV_OBS) $(ODIR)/eigen.o catalog.h
	$(CC) $(CFLAGS) bin_catalog.c $(INCLUDES)  $(CAT_OBJS) $(ODIR)/eigen.o $(SINV_OBS)  \
	-o $(BDIR)/bin_catalog     $(EISLIB) $(LDFLAGS)

$(BDIR)/nsample_catalog: nsample_catalog.c $(CAT_OBJS)  $(SINV_OBS) $(ODIR)/eigen.o catalog.h
	$(CC) $(CFLAGS) nsample_catalog.c $(INCLUDES)  $(CAT_OBJS) $(ODIR)/eigen.o $(SINV_OBS)  \
	-o $(BDIR)/nsample_catalog     $(EISLIB) $(LDFLAGS)

$(BDIR)/solve_stress_one_bin: solve_stress_one_bin.c $(CAT_OBJS)  $(SINV_OBS) $(ODIR)/eigen.o catalog.h
	$(CC) $(CFLAGS) solve_stress_one_bin.c $(INCLUDES)  $(CAT_OBJS) $(ODIR)/eigen.o $(SINV_OBS)  \
	-o $(BDIR)/solve_stress_one_bin     $(EISLIB) $(LDFLAGS)

$(BDIR)/m02mag: m02mag.c $(CAT_OBJS) catalog.h
	$(CC) $(CFLAGS) m02mag.c $(INCLUDES)   $(CAT_OBJS) \
	-o $(BDIR)/m02mag    $(LDFLAGS)

$(BDIR)/m02dcfp: m02dcfp.c $(CAT_OBJS)   catalog.h
	$(CC) $(CFLAGS) m02dcfp.c $(INCLUDES)  $(CAT_OBJS) \
	-o $(BDIR)/m02dcfp      $(LDFLAGS)

$(BDIR)/calc_gr: calc_gr.c  $(CAT_OBJS)   catalog.h
	$(CC) $(CFLAGS) calc_gr.c $(INCLUDES)  $(CAT_OBJS)  \
	-o $(BDIR)/calc_gr      $(LDFLAGS)

$(BDIR)/calc_gr_time: calc_gr_time.c  $(CAT_OBJS)   catalog.h
	$(CC) $(CFLAGS) calc_gr_time.c $(INCLUDES)  $(CAT_OBJS)  \
	-o $(BDIR)/calc_gr_time      $(LDFLAGS)


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

$(BDIR)/test_eigen: test_eigen.c $(CAT_OBJS) $(ODIR)/eigen.o
	$(CC) $(CFLAGS) test_eigen.c $(INCLUDES)  $(CAT_OBJS) \
	-o $(BDIR)/test_eigen  $(ODIR)/eigen.o	$(EISLIB)   $(LDFLAGS)

$(BDIR)/test_kdtree: test_kdtree.c $(ODIR)/geo_kdtree.o $(ODIR)/linalg_misc_geo.o
	$(CC) $(CFLAGS) test_kdtree.c $(INCLUDES) $(ODIR)/geo_kdtree.o $(ODIR)/linalg_misc_geo.o  -o $(BDIR)/test_kdtree $(LDFLAGS)


$(ODIR)/eigen.main.o: eigen_driver.c $(HFILES)
	$(CC) $(CFLAGS) $(DEFINES) -c eigen_driver.c -o $(ODIR)/eigen.main.o 

$(ODIR)/eigen.ov.o: eigen_driver.c $(HFILES)
	$(CC) $(CFLAGS)  $(DEFINES) -c eigen_driver.c -DONLY_VALUES -o $(ODIR)/eigen.ov.o

$(ODIR)/eigen.tds.o: eigen_driver.c $(HFILES)
	$(CC) $(CFLAGS)  $(DEFINES) -c eigen_driver.c -DTHREED_SYMMETRIC -o $(ODIR)/eigen.tds.o

$(ODIR)/eigen.tds.ov.o: eigen_driver.c $(HFILES)
	$(CC) $(CFLAGS) $(DEFINES)  -c eigen_driver.c -DTHREED_SYMMETRIC \
		-DONLY_VALUES -o $(ODIR)/eigen.tds.ov.o


$(ODIR)/%.o: %.c  $(HFILES)
	$(CC) $(CFLAGS) $(DEFINES)  $(INCLUDES) -c $< -o $(ODIR)/$*.o

