LIB_DIR = ..
include $(LIB_DIR)/Makefile_inc
TARGET = $(LIB_DIR)/libem2d.a

OBJS = fio.o particle.o field.o boundary.o mpi_set.o random_gen.o shape_function.o

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update

particle.o : shape_function.o
field.o : boundary.o shape_function.o

clean :
	rm -f $(OBJS) *.mod
