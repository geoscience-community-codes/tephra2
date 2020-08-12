all:    
	cd inversion_src; make -f makefile.linux; cd ..
	cd forward_src; make -f makefile.linux; cd ..
        
.PHONY: clean
clean:  
	cd inversion_src; make -f makefile.linux clean; cd ..
	cd forward_src; make -f makefile.linux clean; cd ..
	rm -vf tephra2_2020 tephra2-inversion_2020
