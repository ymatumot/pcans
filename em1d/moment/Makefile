LIB_DIR = ..
include $(LIB_DIR)/Makefile_inc
TARGET = mom.out
OBJS = boundary.o particle.o const.o main.o mom_calc.o fio.o
.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $< 

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) -o mom.out $(FFLAGS) $(OBJS)

main.o : const.o boundary.o fio.o mom_calc.o particle.o

clean :
	rm -f $(TARGET) $(OBJS) dat/*.dat mom/*.dat psd/*.dat

