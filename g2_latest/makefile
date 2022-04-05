###########################################
# COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999-2000
# FILE NAME: makefile
# CODE CONTROLLER: Jasper Horrell
# DESCRIPTION: 
# Make file for G2 SAR processor. 
#
# VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
# COMMENTS: 
# Add this header info.
#
# VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 2000-03-27
# COMMENTS:
# Add step freq processing module compilation.
#
############################################


CC = gcc
FLAGS = -O2 -Wall 
LIBS = -lm
EXECUTABLES = mocomp rngcom stepf corner azcom sniffdc sunras b2tif imu_unpack iq2mag flt2byte swapend rev4byte pga
RM = rm -f 

all : mocomp rngcom corner azcom sniffdc sunras b2tif imu_unpack iq2mag flt2byte swapend rev4byte pga

g2 : mocomp rngcom stepf corner azcom sniffdc b2tif imu_unpack iq2mag flt2byte swapend pga

swapend : swapend.c
	$(CC) -o swapend swapend.c $(LIBS) $(FLAGS)

iq2mag : iq2mag.c
	$(CC) -o iq2mag iq2mag.c $(LIBS) $(FLAGS) 

imu_unpack : g2imu_unpack.o g2func.o g2parse.o
	$(CC) -o imu_unpack g2imu_unpack.o g2func.o g2parse.o $(LIBS)

g2imu_unpack.o : g2imu_unpack.c g2func.h g2parse.h
	$(CC) -c g2imu_unpack.c $(FLAGS) $(INCLUDE)	

mocomp : g2moc.o g2func.o g2parse.o
	$(CC) -o mocomp g2moc.o g2func.o g2parse.o $(LIBS)

g2moc.o : g2moc.c g2func.h g2moc.h g2parse.h
	$(CC) -c g2moc.c $(FLAGS) $(INCLUDE)	

sniffdc : g2sniffdc.c g2sniffdc.h g2func.h
	$(CC) -o sniffdc g2sniffdc.c $(LIBS) $(FLAGS)

g2func.o : g2func.c g2func.h
	$(CC) -c g2func.c $(FLAGS) $(INCLUDE)

azcom : g2azc.o g2func.o
	$(CC) -o azcom g2azc.o g2func.o $(LIBS)

g2azc.o : g2azc.c g2func.h
	$(CC) -c g2azc.c $(FLAGS) $(INCLUDE) 

rngcom : g2rnc.o g2func.o g2lms.o g2notch.o g2rline.o
	$(CC) -o rngcom g2rnc.o g2func.o g2lms.o g2notch.o g2rline.o $(LIBS)

g2rnc.o : g2rnc.c g2func.h
	$(CC) -c g2rnc.c $(FLAGS) $(INCLUDE) 

g2lms.o : g2lms.c g2func.h
	$(CC) -c g2lms.c $(FLAGS) $(INCLUDE) 

g2notch.o : g2notch.c g2rline.h g2func.h
	$(CC) -c g2notch.c $(FLAGS) $(INCLUDE) 

g2rline.o : g2rline.c g2func.h
	$(CC) -c g2rline.c $(FLAGS) $(INCLUDE) 		

stepf : g2stepf.o g2geth.o sf_procs.o g2func.o
	$(CC) -o stepf g2stepf.o g2geth.o sf_procs.o g2func.o $(LIBS)

g2stepf.o : g2stepf.c g2func.h sf_procs.h g2geth.h
	$(CC) -c g2stepf.c $(FLAGS) $(INCLUDE)

sf_procs.o : sf_procs.c g2func.h
	$(CC) -c sf_procs.c $(FLAGS) $(INCLUDE)

g2geth.o : g2geth.c g2func.h sf_procs.h
	$(CC) -c g2geth.c $(FLAGS) $(INCLUDE)

corner : g2cor.o g2func.o
	$(CC) -o corner g2cor.o g2func.o $(LIBS)

g2cor.o : g2cor.c g2func.h
	$(CC) -c g2cor.c $(FLAGS) $(INCLUDE) 

sunras : g2imgfu.o g2func.o
	$(CC) -o sunras g2imgfu.o g2func.o $(LIBS)

g2imgfu.o : g2imgfu.c g2func.h
	$(CC) -c g2imgfu.c $(FLAGS) $(INCLUDE)  

b2tif : g2b2tif.c g2b2tif.h g2func.h
	$(CC) -o b2tif g2b2tif.c $(LIBS) $(FLAGS)

flt2byte : flt2byte.c
	$(CC) -o flt2byte flt2byte.c $(LIBS) $(FLAGS)

rev4byte : rev4byte.c
	$(CC) -o rev4byte rev4byte.c $(LIBS) $(FLAGS)

g2parse.o : g2parse.c g2func.h g2parse.h
	$(CC) -c g2parse.c $(FLAGS) $(INCLUDE)

pga : pga.o g2func.o
	$(CC) -o pga pga.o g2func.o $(LIBS)

pga.o : pga.c g2func.h
	$(CC) -c pga.c $(FLAGS) $(INCLUDE)

clean:
	$(RM) *.o *.pyc *~ $(EXECUTABLES) core

