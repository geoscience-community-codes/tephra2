all:	
	cd inversion_src; make; cd ..
	cd forward_src; make; cd ..
	
clean:	
	cd inversion_src; make clean; cd ..
	cd forward_src; make clean; cd ..
