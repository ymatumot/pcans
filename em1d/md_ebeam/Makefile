LIB_DIR = ..
include $(LIB_DIR)/Makefile_inc
TARGET = a.out
OBJS = const.o init.o main.o
.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) -L$(LIB_DIR) $(FFLAGS) -I $(LIB_DIR)/common -c $< -lem1d

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) -L$(LIB_DIR) -o $(TARGET) $(FFLAGS) $(OBJS) -lem1d 

$(LIB_DIR)/libem1d.a : 
	cd $(LIB_DIR) ; make

main.o : init.o const.o $(LIB_DIR)/libem1d.a
init.o : const.o $(LIB_DIR)/libem1d.a

clean :
	rm -f $(TARGET) $(OBJS) *.mod dat/*.dat mom/*.dat psd/*.dat

