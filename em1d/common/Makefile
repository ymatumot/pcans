LIB_DIR = ..
include $(LIB_DIR)/Makefile_inc
TARGET = $(LIB_DIR)/libem1d.a

OBJS = fio.o particle.o field.o boundary.o random_gen.o mom_calc.o

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update

field.o : boundary.o
clean :
	rm -f $(OBJS) *.mod
