all:    
#	cd inversion_src; make -f makefile.linux; cd ..
	cd forward_src; make -f makefile.linux; cd ..
        
.PHONY: clean install uninstall
clean:  
#	cd inversion_src; make -f makefile.linux clean; cd ..
	cd forward_src; make -f makefile.linux clean; cd ..

install:
#	cd inversion_src; make -f makefile.linux install; cd ..
	cd forward_src; make -f makefile.linux install; cd ..
	
uninstall:
#	cd inversion_src; make -f makefile.linux uninstall; cd ..
	cd forward_src; make -f makefile.linux uninstall; cd ..

