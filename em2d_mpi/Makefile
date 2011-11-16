LIB_DIR = .
TARGET = $(LIB_DIR)/libem2d.a
all :  $(TARGET)

$(TARGET) : 
	cd common; make

clean :
	cd common; make clean
	cd moment; make clean
	cd md_wave; make clean
	cd md_kh; make clean
	cd md_shock; make clean
	rm -f $(TARGET)

