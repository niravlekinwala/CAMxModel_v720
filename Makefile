#-------------------------------------------------------------------------------
#     This is the Makefile for CAMx v7.20
#
#     Syntax is: 
#     make COMPILER=mycompiler <CONFIG=myapp> <MPI=mpi_option> <NCF=ncf_option> 
#     <IEEE=true>
#
#     The Makefile will generate a CAMx executable program called:
#     CAMx.myapp.mpi_option.<NCF>.mycompiler<.ieee>
#
#     The mandatory "COMPILER" argument defines the compiler to invoke.
#     Generally, to compile CAMx you do not need to use the optional CONFIG
#     argument; if CONFIG is not specified, the default "v7.20" string will be
#     assigned. Depending on if/how the optional MPI argument is set, the
#     "mpi_option" string will be set either to the option you specify or to 
#     "noMPI" (default). The IEEE flag will invoke IEEE-standard math, which
#     will maximize consistency in model results across platforms and compilers.
#     The IEEE flag may cause CAMx to run somewhat slower. The NCF flag
#     will build NetCDF support into the exectuable. This will allow for
#     the option of generating NetCDF files when running the model. The 
#     options are NCF3 (NetCDF version 3), NCF4_C (NetCDF version 4, with 
#     compression) and NCF4_NC (NetCDF version 4, with no compression).
#
#     For example:
#
#     make pgf
#
#     will build a CAMx executable using the Portland Group compiler with OMP
#     and MPI parallelization disabled and using the default include file 
#     "./Includes/camx.prm.v7.20" for configuration parameters. The executable
#     will be named:
#
#     CAMx.v7.20.noMPI.pgf
#
#     As another example:
#
#     make ifortomp CONFIG=mydomain MPI=mpich NCF=NCF4_C IEEE=true
#
#     will build a CAMx executable using the Intel compiler with OMP and MPI
#     parallelization enabled, NetCDF support, IEEE math invoked, and using 
#     the include file "./Includes/camx.prm.mydomain" for the configuration 
#     parameters. The executable will be named:
#
#     CAMx.mydomain.MPICH.NCF4.pgfomp.ieee
#
#  NOTE: RECENT VERSIONS OF THE GFORTRAN COMPILER WILL NOT COMPILE CAMx
#        WITH OMP PARALLELIZATION BECAUSE OF A BLOCK DATA INCOMPATIBILITY IN
#        ISORROPIA.  USE GFORTRAN WITH MPI PARALLELIZATION.  THIS IS NOT AN
#        ISSUE WITH OTHER SUPPORTED COMPILERS.
#-------------------------------------------------------------------------------

.SUFFIXES: .c .f .o .f90 .F90

#-------------------------------------------------------------------------------
# SET PATH TO MPI AND NetCDF INSTALLATIONS
#
#  The NetCDF libraries must be listed in this order:
#    -lnetcdff -lnetcdf
#
MPI_INST = /usr/local/mpich3
NCF_INST = /usr/local/netcdf-4.4.1.1.gfortran

DUM_MPI = ./Dummy_MPI
DUM_NCF = ./Dummy_NCF

#ifndef MPI
   MPI = false
#endif
#ifndef IEEE
   IEEE = false
#endif
#ifndef NCF
   NCF = false
#endif
#
#-------------------------------------------------------------------------------

#ifndef CONFIG
   CONFIG = v7.20
#endif

#-------------------------------------------------------------------------------
# Default CAMx source directories:
#-------------------------------------------------------------------------------

CAMx      = ./CAMx
CMC       = ./CMC
CF_AERO   = ./CF_AERO
CMU_AERO  = ./CMU_AERO
SOAP      = ./SOAP
HG        = ./HG
PIG       = ./PiG
IOBIN     = ./IO_bin
IONCF     = ./IO_NCF
OSAT      = ./OSAT
DDM       = ./DDM
PA        = ./PA
RTRAC     = ./RTRAC
CAMx_INC  = ./Includes
MOD_SRC   = ./Mod_src
MOD_DIR   = ./Modules
CAMX_MPI  = ./MPI

#ifndef COMPILER
   COMPILER = help
#endif

STAR_ECHO = "*****************************************************************"
BLANK_ECHO = "*                                                               *"

#-------------------------------------------------------------------------------
#  MPI specific variables
#-------------------------------------------------------------------------------

ifneq (, $(findstring $(MPI),true mpich MPICH))
    MPI_ECHO = "* MPI will be built in using MPICH1 or MPICH2                   *"
    MPI_INC = $(MPI_INST)/include
    MPI_LIBS = -L$(CAMX_MPI)/util -lparlib \
        -L$(MPI_INST)/lib -lmpich -lpthread
    MPI_STRING = MPICH
else ifneq (, $(findstring $(MPI),mpich3 MPICH3))
    MPI_ECHO = "* MPI will be built in using MPICH3                             *"
    MPI_INC = $(MPI_INST)/include
    MPI_LIBS = -L$(CAMX_MPI)/util -lparlib \
        -L$(MPI_INST)/lib -lmpich -lpthread -lmpl
    MPI_STRING = MPICH3
else ifneq (, $(findstring $(MPI),mvapich MVAPICH))
    MPI_ECHO = "* MPI will be built in using MVAPICH                            *"
    MPI_INC = $(MPI_INST)/include
    MPI_LIBS = -L$(CAMX_MPI)/util -lparlib
    FC = $(MPI_INST)/bin/mpif90
    CC = $(MPI_INST)/bin/mpicc
    MPI_STRING = MVAPICH
else ifneq (, $(findstring $(MPI),openmpi openMPI OPENMPI))
    MPI_ECHO = "* MPI will be built in using OpenMPI                            *"
    MPI_INC = $(MPI_INST)/include
    MPI_LIBS = -L$(CAMX_MPI)/util -lparlib -lnsl -lutil -lm -ldl
    MPI_STRING = openMPI
    FC = $(MPI_INST)/bin/mpif90
    CC = $(MPI_INST)/bin/mpicc
else
    MPI_ECHO = "* MPI will NOT be built in                                      *"
    MPI_INC = $(DUM_MPI)
    MPI_LIBS = 
    MPI_DUM_C = $(DUM_MPI)/mpi_dummy_c.o
    MPI_DUM_F = $(DUM_MPI)/mpi_dummy.o
    MPI_STRING = noMPI
endif

#-------------------------------------------------------------------------------
#  NetCDF specific variables
#-------------------------------------------------------------------------------


ifneq (, $(findstring $(NCF),NCF3 ncf3 Ncf3))
    NCF_ECHO = "* NetCDF will be built in using version 3                       *"
    NCF_INC = $(NCF_INST)/include
    NCF_LIBS = -L$(NCF_INST)/lib -lnetcdf
else ifneq (, $(findstring $(NCF),NCF4_C ncf4_c Ncf4_C Ncf4_c ncf4_C NCF4_c))
    NCF_ECHO = "* NetCDF will be built in using version 4, with compression     *"
    CHUNK = -DCHUNK
    NCF_INC = $(NCF_INST)/include
    NCF_LIBS = -L$(NCF_INST)/lib -lnetcdff -lnetcdf
else ifneq (, $(findstring $(NCF),NCF4_NC ncf4_nc Ncf4_NC Ncf4_Nc ncf4_NC NCF4_nc))
    NCF_ECHO = "* NetCDF will be built in using version 4, with NO compression  *"
    NCF_INC = $(NCF_INST)/include
    NCF_LIBS = -L$(NCF_INST)/lib -lnetcdff -lnetcdf
else
    NCF_ECHO = "* NetCDF will NOT be built in                                   *"
    NCF_INC = $(DUM_NCF)
    NCF_LIBS = 
    NCF_DUM_F = $(DUM_NCF)/ncf_dummy.o
endif

#-------------------------------------------------------------------------------
#  Variables for build
#-------------------------------------------------------------------------------

ifneq (, $(findstring $(IEEE),TRUE true True T))
    IEEE_ECHO = "* The IEEE option will be used                                  *"
    ifneq (, $(findstring $(NCF),NCF3 ncf3 Ncf3))
       TARGT = CAMx.$(CONFIG).$(MPI_STRING).NCF3.ieee.$(COMPILER)
    else ifneq (, $(findstring $(NCF),NCF4_C ncf4_c Ncf4_C Ncf4_c ncf4_C NCF4_c))
       TARGT = CAMx.$(CONFIG).$(MPI_STRING).NCF4.ieee.$(COMPILER)
    else ifneq (, $(findstring $(NCF),NCF4_NC ncf4_nc Ncf4_NC Ncf4_Nc ncf4_NC NCF4_nc))
       TARGT = CAMx.$(CONFIG).$(MPI_STRING).NCF4.ieee.$(COMPILER)
    else
       TARGT = CAMx.$(CONFIG).$(MPI_STRING).ieee.$(COMPILER)
    endif
else
    IEEE_ECHO = "* The IEEE option NOT will be used                              *"
    ifneq (, $(findstring $(NCF),NCF3 ncf3 Ncf3))
        TARGT = CAMx.$(CONFIG).$(MPI_STRING).NCF3.$(COMPILER) 
    else ifneq (, $(findstring $(NCF),NCF4_C ncf4_c Ncf4_C Ncf4_c ncf4_C NCF4_c))
        TARGT = CAMx.$(CONFIG).$(MPI_STRING).NCF4.$(COMPILER) 
    else ifneq (, $(findstring $(NCF),NCF4_NC ncf4_nc Ncf4_NC Ncf4_Nc ncf4_NC NCF4_nc))
       TARGT = CAMx.$(CONFIG).$(MPI_STRING).NCF4.$(COMPILER)
    else
        TARGT = CAMx.$(CONFIG).$(MPI_STRING).$(COMPILER) 
    endif
endif
LIBS     = $(NCF_LIBS) $(MPI_LIBS)
INCLUDES = -I$(CAMx_INC) -I$(MPI_INC) -I$(NCF_INC)

#-------------------------------------------------------------------------------
#  Compiler specific variables
#-------------------------------------------------------------------------------

ifneq (, $(findstring $(COMPILER),pgf pgfomp))
    FC = pgf90
    CC = gcc
    ifneq (, $(findstring $(MPI_STRING),mpich MPICH mpich3 MPICH3 openmpi openMPI mvapich OPENMPI MVAPICH))
       FC = $(MPI_INST)/bin/mpif90
       CC = $(MPI_INST)/bin/mpicc
    endif
    FLGS = $(INCLUDES) -mcmodel=medium -O2 -Mdalign -Mextend -Mnoframe -byteswapio -Bstatic_pgi
    ifeq ($(COMPILER),pgfomp)
       FLGS_OMP = -mp
    endif
    ifneq (, $(findstring $(IEEE),TRUE true True T))
       FLGS_IEEE = -Kieee
    endif
    MODULES=-I${MOD_DIR} -module ${MOD_DIR}
else ifneq (, $(findstring $(COMPILER),ifort ifortomp))
    FC = ifort
    CC = gcc
    ifneq (, $(findstring $(MPI_STRING),openmpi openMPI mvapich OPENMPI MVAPICH))
       FC = $(MPI_INST)/bin/mpif90
       CC = $(MPI_INST)/bin/mpicc
    endif
    FLGS = $(INCLUDES) -mcmodel=medium -O2 -align dcommons -extend_source -convert big_endian
    MODULES=-I${MOD_DIR} -module ${MOD_DIR}
    ifeq ($(COMPILER),ifortomp)
       FLGS_OMP = -openmp
    endif
    ifneq (, $(findstring $(IEEE),TRUE true True T))
       FLGS_IEEE = -mieee-fp
    endif
else ifneq (, $(findstring $(COMPILER),absoft absoftomp))
    FC = af90
    CC = cc
    FLGS=$(INCLUDES) -V -O2 -m64 -YEXT_NAMES=LCS -YEXT_SFX=_ -s -YCFRL=1 -W132 -lU77 -bigendian -L/Applications/Absoft/lib64
    MODULES=-I${MOD_DIR}
    ifeq ($(COMPILER),absoftomp)
       FLGS_OMP = -openmp
    endif
else ifneq (, $(findstring $(COMPILER),gfortran gfortranomp))
    FC = gfortran
    CC = gcc
    FLGS=$(INCLUDES) -O2 -mcmodel=medium -J${MOD_DIR} -fno-align-commons -fconvert=big-endian -frecord-marker=4 -ffixed-line-length-0
    MODULES=-I${MOD_DIR}
    ifneq (, $(findstring $(MPI_STRING),mpich MPICH mpich3 MPICH3 openmpi openMPI mvapich OPENMPI MVAPICH))
       FC = $(MPI_INST)/bin/mpif90
       CC = $(MPI_INST)/bin/mpicc
       MPI_LIBS = -L$(CAMX_MPI)/util -lparlib
    endif
    ifeq ($(COMPILER),gfortranomp)
       FLGS_OMP = -fopenmp
    endif
else ifneq (, $(findstring $(COMPILER),oracle oracleomp))
    FC = f90
    CC = gcc
    ifneq (, $(findstring $(MPI_STRING),openmpi openMPI mvapich OPENMPI MVAPICH))
       FC = $(MPI_INST)/bin/mpif90
       CC = $(MPI_INST)/bin/mpicc
    endif
    FLGS=$(INCLUDES) -O3 -aligncommon -e -xfilebyteorder=big4:%all
    MODULES=-M${MOD_DIR}
    ifneq (, $(findstring $(MPI_STRING),mpich MPICH mpich3 MPICH3 openmpi openMPI mvapich OPENMPI MVAPICH))
       FC = $(MPI_INST)/bin/mpif90
       CC = $(MPI_INST)/bin/mpicc
       MPI_LIBS = -L$(CAMX_MPI)/util -lparlib
    endif
    ifeq ($(COMPILER),oracleomp)
       FLGS_OMP = -openmp
    endif
endif

#-------------------------------------------------------------------------------
#  Main directive
#-------------------------------------------------------------------------------

all: comp_$(COMPILER)

#-------------------------------------------------------------------------------
#  Directives for error trapping
#-------------------------------------------------------------------------------

pgf: comp_help
pgfomp: comp_help
ifort: comp_help
ifortomp: comp_help

pg_linux: comp_help
pg_linuxomp: comp_help
i_linux: comp_help
i_linuxomp: comp_help
absoft: comp_help
absoftomp: comp_help
comp_pg_linux: comp_help
comp_pg_linuxomp: comp_help
comp_i_linux: comp_help
comp_i_linuxomp: comp_help
comp_absoft: comp_help
comp_absoftomp: comp_help

help:
	@echo '----------------------------------------------------------------'
	@echo ' '
	@echo 'This is the Makefile for CAMx v7.20'
	@echo ' '
	@echo 'Syntax is:'
	@echo 'make COMPILER=mycompiler <CONFIG=myapp> <MPI=mpi_option> <NCF=ncf_option> <IEEE=true>'
	@echo ' '
	@echo 'Acceptable entries for COMPILER are:'
	@echo '    pgf         -- PGF90 for Linux'
	@echo '    pgfomp      -- PGF90 w/OMP for Linux'
	@echo '    ifort       -- IFORT for Linux'
	@echo '    ifortomp    -- IFORT w/OMP for Linux'
	@echo '    gfortran    -- GFortran for Linux'
	@echo '    gfortranomp -- GFortran w/OMP for Linux'
	@echo '    oracle      -- Oracle Sun Studio for Linux'
	@echo '    oracleomp   -- Oracle Sun Studio w/OMP for Linux'
	@echo '    absoft      -- ABSoft for Mac OSX'
	@echo '    absoftomp   -- ABSoft w/OMP for Mac OSX'
	@echo ' '
	@echo 'CONFIG=myapp points the make utility at the specific'
	@echo 'CAMx parameters file called "./Includes/camx.prm.myapp"'
	@echo 'which you have created for your specific CAMx application.'
	@echo 'If CONFIG is not set on the command line, the default'
	@echo 'CAMx parameters file will be used: "./Includes/camx.prm.v7.20"'
	@echo ' '
	@echo 'MPI=mpi_option will compile CAMx to include MPI so that the'
	@echo 'model can be run in parallel on distributed memory systems. MPI'
	@echo 'can be run in tandem with OMP parallelization on shared memory'
	@echo 'systems. Acceptable entries for MPI are:'
	@echo '    mpich  (for MPICH1/2)'
	@echo '    mpich3 '
	@echo '    mvapich'
	@echo '    openmpi (PGF and IFORT only)'
	@echo 'Your system must have the specified MPI support installed and '
	@echo 'you may need to change the MPI_INST variable in the Makefile '
	@echo 'to the path of the MPI installation.'
	@echo ' '
	@echo 'The NCF flag will build NetCDF support into the model. You'
	@echo 'will need to ensure that the NCF_INST makefile variable is'
	@echo 'properly defined. With NetCDF support you will have the '
	@echo 'option at runtime to generate NetCDF output files.'
	@echo 'Acceptable entries for NCF are:'
	@echo '    NCF3    - NetCDF support with NetCDF libraries from version 3.'
	@echo '    NCF4_C  - NetCDF support with NetCDF libraries from version 4 that'
	@echo '              include compression (HDF5 and Zlib)'
	@echo '    NCF4_NC - NetCDF support with NetCDF libraries from version 4 that'
	@echo '              do not include compression'
	@echo '    False   - no NETCD support (default)'
	@echo ' '
	@echo 'The IEEE flag will invoke IEEE-standard math, which will'
	@echo 'maximize consistency in model results across platforms and'
	@echo 'compilers. The IEEE flag may cause CAMx to run somewhat slower.'
	@echo 'The IEEE option is available only for PGF and IFORT compilers.'
	@echo ' '
	@echo '----------------------------------------------------------------'

comp_help:
	@echo '----------------------------------------------------------------'
	@echo ' '
	@echo 'This is the Makefile for CAMx v7.20'
	@echo ' '
	@echo 'Syntax is:'
	@echo 'make COMPILER=mycompiler <CONFIG=myapp> <MPI=mpi_option> <NCF=ncf_option> <IEEE=true>'
	@echo ' '
	@echo 'Acceptable entries for COMPILER are:'
	@echo '    pgf         -- PGF90 for Linux'
	@echo '    pgfomp      -- PGF90 w/OMP for Linux'
	@echo '    ifort       -- IFORT for Linux'
	@echo '    ifortomp    -- IFORT w/OMP for Linux'
	@echo '    gfortran    -- GFortran for Linux'
	@echo '    gfortranomp -- GFortran w/OMP for Linux'
	@echo '    oracle      -- Oracle Sun Studio for Linux'
	@echo '    oracleomp   -- Oracle Sun Studio w/OMP for Linux'
	@echo '    absoft      -- ABSoft for Mac OSX'
	@echo '    absoftomp   -- ABSoft w/OMP for Mac OSX'
	@echo ' '
	@echo 'CONFIG=myapp points the make utility at the specific'
	@echo 'CAMx parameters file called "./Includes/camx.prm.myapp"'
	@echo 'which you have created for your specific CAMx application.'
	@echo 'If CONFIG is not set on the command line, the default'
	@echo 'CAMx parameters file will be used: "./Includes/camx.prm.v7.20"'
	@echo ' '
	@echo 'MPI=mpi_option will compile CAMx to include MPI so that the'
	@echo 'model can be run in parallel on distributed memory systems. MPI'
	@echo 'can be run in tandem with OMP parallelization on shared memory'
	@echo 'systems. Acceptable entries for MPI are:'
	@echo '    mpich  (for MPICH1/2)'
	@echo '    mpich3 '
	@echo '    mvapich'
	@echo '    openmpi (PGF and IFORT only)'
	@echo 'Your system must have the specified MPI support installed and '
	@echo 'you may need to change the MPI_INST variable in the Makefile '
	@echo 'to the path of the MPI installation.'
	@echo ' '
	@echo 'The IEEE flag will invoke IEEE-standard math, which will'
	@echo 'maximize consistency in model results across platforms and'
	@echo 'compilers. The IEEE flag may cause CAMx to run somewhat slower.'
	@echo 'The IEEE option is available only for PGF and IFORT compilers.'
	@echo ' '
	@echo 'The NCF flag will build NetCDF support into the model. You'
	@echo 'will need to ensure that the NCF_INST makefile variable is'
	@echo 'properly defined. With NetCDF support you will have the '
	@echo 'option at runtime to generate NetCDF output files.'
	@echo 'Acceptable entries for NCF are:'
	@echo '    NCF3    - NetCDF support with NetCDF libraries from version 3.'
	@echo '    NCF4_C  - NetCDF support with NetCDF libraries from version 4 that'
	@echo '              include compression (HDF5 and Zlib)'
	@echo '    NCF4_NC - NetCDF support with NetCDF libraries from version 4 that'
	@echo '              do not include compression'
	@echo '    False   - no NETCD support (default)'
	@echo ' '
	@echo '----------------------------------------------------------------'

#-------------------------------------------------------------------------------
#  Configuration specific directives
#-------------------------------------------------------------------------------

comp_pgf:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_pgfomp:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_ifort:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_ifortomp:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_absoft:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_absoftomp:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_gfortran:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_gfortranomp:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_oracle:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

comp_oracleomp:
	@echo $(STAR_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(MPI_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(NCF_ECHO)
	@echo $(BLANK_ECHO)
	@echo $(IEEE_ECHO)
	@echo $(BLANK_ECHO)
	@echo "*      Executable is: "$(TARGT)
	@echo $(BLANK_ECHO)
	@echo $(STAR_ECHO)
	@sleep 10
	@rm -f $(CAMx_INC)/camx.prm
	@ln -s camx.prm.$(CONFIG) $(CAMx_INC)/camx.prm
	make modules
	make model

clean:	
	rm -f $(OBJCTS) $(MOD_OBJCTS) $(DUM_MPI)/*.o $(DUM_NCF)/*.o *.mod $(MOD_DIR)/*.mod $(CAMX_MPI)/*.o *.a
	cd MPI/util; make clean


MOD_OBJCTS = \
$(MOD_SRC)/camx_includes.o \
$(MOD_SRC)/grid_dims_mod.o \
$(MOD_SRC)/node_mod.o \
$(MOD_SRC)/o3colmap.o \
$(MOD_SRC)/bndary.o \
$(MOD_SRC)/camxcom.o \
$(MOD_SRC)/camxfld.o \
$(MOD_SRC)/chmstry.o \
$(MOD_SRC)/filunit.o \
$(MOD_SRC)/grid.o \
$(MOD_SRC)/grid_nodes.o \
$(MOD_SRC)/master_mod.o \
$(MOD_SRC)/ncf_iomod.o \
$(MOD_SRC)/pigsty.o \
$(MOD_SRC)/ptemiss.o \
$(MOD_SRC)/procan.o \
$(MOD_SRC)/rtracchm.o \
$(MOD_SRC)/rtcmcchm.o \
$(MOD_SRC)/tracer.o

OBJCTS = \
$(CAMx)/CAMx.o \
$(CAMx)/aerochem_cf.o \
$(CAMx)/aerochem_cmu.o \
$(CAMx)/aeroset.o \
$(CAMx)/aggdep.o \
$(CAMx)/aggr00.o \
$(CAMx)/aggreg.o \
$(CAMx)/banner.o \
$(CAMx)/o3colprep.o \
$(CAMx)/o3col_upd.o \
$(CAMx)/average.o \
$(CAMx)/avgall.o \
$(CAMx)/bc1grd.o \
$(CAMx)/bcmodfy.o \
$(CAMx)/bndry_updt.o \
$(CAMx)/caldate.o \
$(CAMx)/camxerr.o \
$(CAMx)/chem10.o \
$(CAMx)/chemdriv.o \
$(CAMx)/chemrxn.o \
$(CAMx)/chmdat.o \
$(CAMx)/chrtime.o \
$(CAMx)/cigdrive.o \
$(CAMx)/cpivot.o \
$(CAMx)/convtrans.o \
$(CAMx)/cvtdate.o \
$(CAMx)/cvtwind.o \
$(CAMx)/depsmry.o \
$(CAMx)/diffus.o \
$(CAMx)/dlsode.o \
$(CAMx)/drvtuv.o \
$(CAMx)/drydep.o \
$(CAMx)/ebisolv.o \
$(CAMx)/emassign.o \
$(CAMx)/emiss.o \
$(CAMx)/emistrns.o \
$(CAMx)/emiss_updt.o \
$(CAMx)/exptbl.o \
$(CAMx)/fgavrg.o \
$(CAMx)/finwind.o \
$(CAMx)/getalbedo.o \
$(CAMx)/findtrop.o \
$(CAMx)/fndkey.o \
$(CAMx)/getdelt.o \
$(CAMx)/getime.o \
$(CAMx)/getunit.o \
$(CAMx)/getznth.o \
$(CAMx)/grdprep.o \
$(CAMx)/hadvbot.o \
$(CAMx)/hadvppm.o \
$(CAMx)/henryfnc.o \
$(CAMx)/iassgn2d.o \
$(CAMx)/iniptr.o \
$(CAMx)/initnml.o \
$(CAMx)/interp2d.o \
$(CAMx)/intrpcnc.o \
$(CAMx)/intrpdat.o \
$(CAMx)/ixemis.o \
$(CAMx)/istrln.o \
$(CAMx)/jstlft.o \
$(CAMx)/juldate.o \
$(CAMx)/khetero.o \
$(CAMx)/khorz.o \
$(CAMx)/kphoto.o \
$(CAMx)/ktherm.o \
$(CAMx)/lcpgeo.o \
$(CAMx)/linpac.o \
$(CAMx)/lscall.o \
$(CAMx)/lsode.o \
$(CAMx)/luassgn.o \
$(CAMx)/massum.o \
$(CAMx)/micromet.o \
$(CAMx)/mrcgeo.o \
$(CAMx)/nesting.o \
$(CAMx)/nstprep.o \
$(CAMx)/parntchd.o \
$(CAMx)/plumeris.o \
$(CAMx)/pntprep.o \
$(CAMx)/pnt_update.o \
$(CAMx)/pspgeo.o \
$(CAMx)/raddrivr.o \
$(CAMx)/rassgn2d.o \
$(CAMx)/rassgn3d.o \
$(CAMx)/rassgn4d.o \
$(CAMx)/reado3col.o \
$(CAMx)/readchm.o \
$(CAMx)/readnml.o \
$(CAMx)/readpht.o \
$(CAMx)/rpsgeo.o \
$(CAMx)/scavrat.o \
$(CAMx)/setbc1d.o \
$(CAMx)/setbc.o \
$(CAMx)/sim_init.o \
$(CAMx)/srfmod.o \
$(CAMx)/srfruf.o \
$(CAMx)/startup.o \
$(CAMx)/strato3.o \
$(CAMx)/timestep.o \
$(CAMx)/timrates.o \
$(CAMx)/toupper.o \
$(CAMx)/tridiag.o \
$(CAMx)/tstep_init.o \
$(CAMx)/tuv.o \
$(CAMx)/updtmet.o \
$(CAMx)/uptime.o \
$(CAMx)/utmgeo.o \
$(CAMx)/vadvppm.o \
$(CAMx)/vd_aer.o \
$(CAMx)/vd_aer_zhang.o \
$(CAMx)/vd_gas.o \
$(CAMx)/vd_gas_zhang.o \
$(CAMx)/vdiffacmx.o \
$(CAMx)/vdiffk.o \
$(CAMx)/ve_gas.o \
$(CAMx)/vnmshcal.o \
$(CAMx)/vrtslv.o \
$(CAMx)/wetdep.o \
$(CAMx)/wrtmass.o \
$(CAMx)/xerf.o \
$(CAMx)/xyadvec.o \
$(CAMx)/zadvec.o \
$(CAMx)/zeros.o \
$(CAMx)/zrates.o \
$(IOBIN)/areaprep.o \
$(IOBIN)/bndprep.o \
$(IOBIN)/clciwt.o \
$(IOBIN)/cncprep.o \
$(IOBIN)/depprep.o \
$(IOBIN)/get_nptsrc.o \
$(IOBIN)/getdepth.o \
$(IOBIN)/openfils.o \
$(IOBIN)/metinit.o \
$(IOBIN)/rdfgcon.o \
$(IOBIN)/rdmethdr.o \
$(IOBIN)/rdpthdr.o \
$(IOBIN)/rdsrfmod.o \
$(IOBIN)/rdsumbc.o \
$(IOBIN)/readar.o \
$(IOBIN)/readbnd.o \
$(IOBIN)/readcnc.o \
$(IOBIN)/readinp.o \
$(IOBIN)/readpt.o \
$(IOBIN)/readtop.o \
$(IOBIN)/srfprep.o \
$(IOBIN)/topprep.o \
$(IOBIN)/wrfgcon.o \
$(IOBIN)/wrtcon.o \
$(IOBIN)/wrtdep.o \
$(IOBIN)/wrtsrf.o \
$(IONCF)/ncf_areaprep.o \
$(IONCF)/ncf_bndprep.o \
$(IONCF)/ncf_chkfile.o \
$(IONCF)/ncf_chk_griddef.o \
$(IONCF)/ncf_chk_tstep.o \
$(IONCF)/ncf_clciwt.o \
$(IONCF)/ncf_closefiles.o \
$(IONCF)/ncf_cncprep.o \
$(IONCF)/ncf_createfile.o \
$(IONCF)/ncf_enddef_file.o \
$(IONCF)/ncf_getdepth.o \
$(IONCF)/ncf_get_tstep.o \
$(IONCF)/ncf_get_nlayers.o \
$(IONCF)/ncf_luseprep.o \
$(IONCF)/ncf_metinit_2d.o \
$(IONCF)/ncf_metinit_3d.o \
$(IONCF)/ncf_metinit_kv.o \
$(IONCF)/ncf_metprep.o \
$(IONCF)/ncf_rdsumbc.o \
$(IONCF)/ncf_readar.o \
$(IONCF)/ncf_readbnd.o \
$(IONCF)/ncf_readcnc.o \
$(IONCF)/ncf_readinp_2d.o \
$(IONCF)/ncf_readinp_3d.o \
$(IONCF)/ncf_readinp_cld.o \
$(IONCF)/ncf_readinp_kv.o \
$(IONCF)/ncf_readinp_lu.o \
$(IONCF)/ncf_readpt.o \
$(IONCF)/ncf_readtop.o \
$(IONCF)/ncf_rdpthdr.o \
$(IONCF)/ncf_rdstacks.o \
$(IONCF)/ncf_set_cache.o \
$(IONCF)/ncf_set_compress_flag.o \
$(IONCF)/ncf_set_emiss_mapping.o \
$(IONCF)/ncf_set_global.o \
$(IONCF)/ncf_set_global_smp.o \
$(IONCF)/ncf_set_specatt_avrg.o \
$(IONCF)/ncf_set_specatt_cpa.o \
$(IONCF)/ncf_set_specatt_ddm.o \
$(IONCF)/ncf_set_specatt_dep.o \
$(IONCF)/ncf_set_specatt_rtrac.o \
$(IONCF)/ncf_set_specatt_rtrac_dep.o \
$(IONCF)/ncf_set_specatt_sadep.o \
$(IONCF)/ncf_set_specatt_sa.o \
$(IONCF)/ncf_set_specatt_srfmod.o \
$(IONCF)/ncf_set_species_mapping.o \
$(IONCF)/ncf_set_tstep.o \
$(IONCF)/ncf_set_vars_base.o \
$(IONCF)/ncf_topprep.o \
$(IONCF)/ncf_wrt_data_dep.o \
$(IONCF)/ncf_wrt_data_grid.o \
$(IONCF)/ncf_wrt_data_gridsmp.o \
$(IONCF)/ncf_wrt_data_species.o \
$(IONCF)/ncf_wrt_data_srfmod.o \
$(IONCF)/ncf_wrt_data_tstep.o \
$(IONCF)/ncf_wrt_dim.o \
$(IONCF)/ncf_wrt_global_ddm.o \
$(IONCF)/ncf_wrt_global.o \
$(IONCF)/ncf_wrt_global_sa.o \
$(IONCF)/ncf_wrt_vars_base.o \
$(IONCF)/ncf_wrt_vars_species.o \
$(IONCF)/ncf_wrt_topo.o \
$(CF_AERO)/hlconst.o \
$(CF_AERO)/hlindex.o \
$(CF_AERO)/isocom_v1.7.o \
$(CF_AERO)/isofwd_v1.7.o \
$(CF_AERO)/isorev_v1.7.o \
$(CF_AERO)/raqchem.o \
$(CF_AERO)/eqsam4clim.o \
$(CMC)/ddmrate1.o \
$(CMC)/ddmrate3.o \
$(CMC)/ddmrate4.o \
$(CMC)/ddmrate5.o \
$(CMC)/ddmrate6.o \
$(CMC)/ddmrate7.o \
$(CMC)/ebirate1.o \
$(CMC)/ebirate3.o \
$(CMC)/ebirate4.o \
$(CMC)/ebirate5.o \
$(CMC)/ebirate6.o \
$(CMC)/ebirate7.o \
$(CMC)/ebirxn1.o \
$(CMC)/ebirxn3.o \
$(CMC)/ebirxn4.o \
$(CMC)/ebirxn5.o \
$(CMC)/ebirxn6.o \
$(CMC)/ebirxn7.o \
$(CMC)/hddmjac1.o \
$(CMC)/hddmjac3.o \
$(CMC)/hddmjac4.o \
$(CMC)/hddmjac5.o \
$(CMC)/hddmjac6.o \
$(CMC)/hddmjac7.o \
$(CMC)/hr_hox1.o \
$(CMC)/hr_hox3.o \
$(CMC)/hr_hox4.o \
$(CMC)/hr_hox5.o \
$(CMC)/hr_hox6.o \
$(CMC)/hr_hox7.o \
$(CMC)/hr_nox1.o \
$(CMC)/hr_nox3.o \
$(CMC)/hr_nox4.o \
$(CMC)/hr_nox5.o \
$(CMC)/hr_nox6.o \
$(CMC)/hr_nox7.o \
$(CMC)/hr_nxy1.o \
$(CMC)/hr_nxy3.o \
$(CMC)/hr_nxy4.o \
$(CMC)/hr_nxy5.o \
$(CMC)/hr_nxy6.o \
$(CMC)/hr_nxy7.o \
$(CMC)/hr_pan1.o \
$(CMC)/hr_pan3.o \
$(CMC)/hr_pan4.o \
$(CMC)/hr_pan5.o \
$(CMC)/hr_pan6.o \
$(CMC)/hr_pan7.o \
$(CMC)/hr_xo1.o \
$(CMC)/hr_xo3.o \
$(CMC)/hr_xo4.o \
$(CMC)/hr_xo7.o \
$(CMC)/hr_xo_dum.o \
$(CMC)/lsrate_pig.o \
$(CMC)/lsrate1.o \
$(CMC)/lsrate3.o \
$(CMC)/lsrate4.o \
$(CMC)/lsrate5.o \
$(CMC)/lsrate6.o \
$(CMC)/lsrate7.o \
$(CMC)/lsrxn1.o \
$(CMC)/lsrxn3.o \
$(CMC)/lsrxn4.o \
$(CMC)/lsrxn5.o \
$(CMC)/lsrxn6.o \
$(CMC)/lsrxn7.o \
$(CMC)/lsrxn_pig.o \
$(CMU_AERO)/areag.o \
$(CMU_AERO)/addit.o \
$(CMU_AERO)/aerchem.o \
$(CMU_AERO)/aero_err.o \
$(CMU_AERO)/aqchem.o \
$(CMU_AERO)/aqdist.o \
$(CMU_AERO)/aqfex1.o \
$(CMU_AERO)/aqfex2.o \
$(CMU_AERO)/aqintegr1.o \
$(CMU_AERO)/aqintegr2.o \
$(CMU_AERO)/aqjex.o \
$(CMU_AERO)/aqoperator1.o \
$(CMU_AERO)/aqrates1.o \
$(CMU_AERO)/aqrates2.o \
$(CMU_AERO)/baddit.o \
$(CMU_AERO)/block.o \
$(CMU_AERO)/bmass.o \
$(CMU_AERO)/coagul.o \
$(CMU_AERO)/constants.o \
$(CMU_AERO)/decisions.o \
$(CMU_AERO)/diameter.o \
$(CMU_AERO)/diff.o \
$(CMU_AERO)/differ.o \
$(CMU_AERO)/diffund.o \
$(CMU_AERO)/dropinit.o \
$(CMU_AERO)/electro.o \
$(CMU_AERO)/eqpart.o \
$(CMU_AERO)/eqparto.o \
$(CMU_AERO)/equaer.o \
$(CMU_AERO)/fullequil.o \
$(CMU_AERO)/hybrd.o \
$(CMU_AERO)/hybrid.o \
$(CMU_AERO)/linint.o \
$(CMU_AERO)/madm.o \
$(CMU_AERO)/mass.o \
$(CMU_AERO)/negchk.o \
$(CMU_AERO)/newdist.o \
$(CMU_AERO)/nucl.o \
$(CMU_AERO)/qsaturation.o \
$(CMU_AERO)/react.o \
$(CMU_AERO)/state.o \
$(CMU_AERO)/steady.o \
$(CMU_AERO)/step.o \
$(CMU_AERO)/svode.o \
$(CMU_AERO)/values.o \
$(CMU_AERO)/vsrm.o \
$(CMU_AERO)/slsode.o \
$(DDM)/adjbcddm.o \
$(DDM)/adjddmc0.o \
$(DDM)/avgrcpddm.o \
$(DDM)/bottddm.o \
$(DDM)/cvticddm.o \
$(DDM)/ddmebi.o \
$(DDM)/ddmeqnslv.o \
$(DDM)/ddmisrpia.o \
$(DDM)/ddmradm.o \
$(DDM)/ddmsoap.o \
$(DDM)/filspddm.o \
$(DDM)/hddmadj.o \
$(DDM)/hddmchem.o \
$(DDM)/hdrrcpddm.o \
$(DDM)/loaddm.o \
$(DDM)/loaddmpm.o \
$(DDM)/ncf_rdicddm.o \
$(DDM)/ncf_rdbcddm.o \
$(DDM)/ncf_rdtcddm.o \
$(DDM)/ncf_rdpthdr_ddm.o \
$(DDM)/ncf_rdptddm.o \
$(DDM)/ncf_emprepddm.o \
$(DDM)/ncf_rdarddm.o \
$(DDM)/rdarddm.o \
$(DDM)/rdbcddm.o \
$(DDM)/rdtcddm.o \
$(DDM)/rdicddm.o \
$(DDM)/rdoptddm.o \
$(DDM)/rdptddm.o \
$(DDM)/specddm.o \
$(DDM)/startddm.o \
$(OSAT)/chkpartial.o \
$(OSAT)/addrcp.o \
$(OSAT)/adjstc0.o \
$(OSAT)/avgrcp.o \
$(OSAT)/avgwal.o \
$(OSAT)/caliball.o \
$(OSAT)/clcbwt.o \
$(OSAT)/clcewt.o \
$(OSAT)/clrbdysa.o \
$(OSAT)/cyctpnsa.o \
$(OSAT)/emprepsa.o \
$(OSAT)/filaqsa.o \
$(OSAT)/filbdysa.o \
$(OSAT)/filptopsa.o \
$(OSAT)/filvdsa.o \
$(OSAT)/get_nptsrc_sa.o \
$(OSAT)/hdrrcp.o \
$(OSAT)/hdrdepsa.o \
$(OSAT)/hdrwsa.o \
$(OSAT)/initsa.o \
$(OSAT)/loadinstsa.o \
$(OSAT)/ncf_emprepsa.o \
$(OSAT)/ncf_rdargrp.o \
$(OSAT)/ncf_rdstacks_sa.o \
$(OSAT)/ncf_rdpthdr_sa.o \
$(OSAT)/ncf_rdptgrp.o \
$(OSAT)/ncf_rd_sa_bcfile.o \
$(OSAT)/ncf_rd_sa_icfile.o \
$(OSAT)/ncf_rd_sa_tcfile.o \
$(OSAT)/ncf_sumgrps.o \
$(OSAT)/ncf_sa_bcprep.o \
$(OSAT)/ncf_sa_icprep.o \
$(OSAT)/ncf_sa_tcprep.o \
$(OSAT)/osatsa.o \
$(OSAT)/o3prdsa.o \
$(OSAT)/psatsa.o \
$(OSAT)/pigdrysa.o \
$(OSAT)/pigdumpsa.o \
$(OSAT)/pigwetsa.o \
$(OSAT)/rdargrp.o \
$(OSAT)/rdfgsa.o \
$(OSAT)/rdinstsa.o \
$(OSAT)/rdoptsa.o \
$(OSAT)/rdpartial_map.o \
$(OSAT)/rdpthdr_sa.o \
$(OSAT)/rdptgrp.o \
$(OSAT)/readarsa.o \
$(OSAT)/readptsa.o \
$(OSAT)/recalib.o \
$(OSAT)/recalib2.o \
$(OSAT)/repartsa.o \
$(OSAT)/rercp.o \
$(OSAT)/sa_bcprep.o \
$(OSAT)/sa_icprep.o \
$(OSAT)/sa_tcprep.o \
$(OSAT)/spcsprcsa.o \
$(OSAT)/spccb05sa.o \
$(OSAT)/spccb6r4sa.o \
$(OSAT)/spccb7sa.o \
$(OSAT)/rd_sa_bcfile.o \
$(OSAT)/rd_sa_icfile.o \
$(OSAT)/rd_sa_tcfile.o \
$(OSAT)/resmap.o \
$(OSAT)/specsa.o \
$(OSAT)/stabsa.o \
$(OSAT)/startsa.o \
$(OSAT)/sum1grd.o \
$(OSAT)/sum1pnt.o \
$(OSAT)/sumgrps.o \
$(OSAT)/sumicwt.o \
$(OSAT)/sumwt4.o \
$(OSAT)/timadv.o \
$(OSAT)/wconsa.o \
$(OSAT)/wfconsa.o \
$(OSAT)/wrtdepsa.o \
$(OSAT)/wsfcsa.o \
$(OSAT)/xfluxsa.o \
$(OSAT)/yfluxsa.o \
$(PA)/cpamech1.o \
$(PA)/cpamech3.o \
$(PA)/cpamech4.o \
$(PA)/cpamech5.o \
$(PA)/cpamech6.o \
$(PA)/cpamech7.o \
$(PA)/cparad.o \
$(PA)/initipr.o \
$(PA)/paconv.o \
$(PA)/pagrids.o \
$(PA)/pasetup.o \
$(PA)/pazero.o \
$(PA)/rdoptpa.o \
$(PA)/wrtcpa.o \
$(PA)/wrtipr.o \
$(PA)/wrtiprhdr.o \
$(PA)/wrtirr.o \
$(PA)/wrtirrhdr.o \
$(PIG)/avepig.o \
$(PIG)/irondry.o \
$(PIG)/ironchem.o \
$(PIG)/ironwet.o \
$(PIG)/noxchem.o \
$(PIG)/pigcoord.o \
$(PIG)/pigdiag.o \
$(PIG)/pigdrive.o \
$(PIG)/pigevol.o \
$(PIG)/piggrow.o \
$(PIG)/piginit.o \
$(PIG)/pigmscl.o \
$(PIG)/pigprep.o \
$(PIG)/pigwalk.o \
$(PIG)/pigsampl.o \
$(PIG)/siggrow.o \
$(PIG)/smpl1puf.o \
$(PIG)/smpprep.o \
$(PIG)/virtdump.o \
$(PIG)/walk1puf.o \
$(PIG)/walk1pufz.o \
$(PIG)/wrtpig.o \
$(PIG)/wrtsmp.o \
$(RTRAC)/chemrt.o \
$(RTRAC)/cvticrt.o \
$(RTRAC)/deqmcmc.o \
$(RTRAC)/djaccmc.o \
$(RTRAC)/dlsdrv.o \
$(RTRAC)/dprepslo.o \
$(RTRAC)/dratecmc.o \
$(RTRAC)/drydeprt.o \
$(RTRAC)/dslocmc.o \
$(RTRAC)/empreprt.o \
$(RTRAC)/gencmc.o \
$(RTRAC)/get_nptsrc_rt.o \
$(RTRAC)/hdrcprt.o \
$(RTRAC)/hdrwsrf.o \
$(RTRAC)/idxspc.o \
$(RTRAC)/is1alpha.o \
$(RTRAC)/is1num.o \
$(RTRAC)/krtc.o \
$(RTRAC)/ksci.o \
$(RTRAC)/ncf_empreprt.o \
$(RTRAC)/ncf_rdarrt.o \
$(RTRAC)/ncf_rdbcrt.o \
$(RTRAC)/ncf_rdicrt.o \
$(RTRAC)/ncf_rdptrt.o \
$(RTRAC)/partitionrt.o \
$(RTRAC)/pigchmrt.o \
$(RTRAC)/rbjaccmc.o \
$(RTRAC)/rbkdrv.o \
$(RTRAC)/rbratcmc.o \
$(RTRAC)/rdarrt.o \
$(RTRAC)/rdbcrt.o \
$(RTRAC)/rdchmrt.o \
$(RTRAC)/rdicrt.o \
$(RTRAC)/rdozmch.o \
$(RTRAC)/rdptrt.o \
$(RTRAC)/rdrcprt.o \
$(RTRAC)/rdschm.o \
$(RTRAC)/rdsrfrt.o \
$(RTRAC)/reemisrt.o \
$(RTRAC)/rodas.o \
$(RTRAC)/seqmcmc.o \
$(RTRAC)/sjaccmc.o \
$(RTRAC)/slsdrv.o \
$(RTRAC)/sprepslo.o \
$(RTRAC)/sratecmc.o \
$(RTRAC)/srfmodrt.o \
$(RTRAC)/sslocmc.o \
$(RTRAC)/startrt.o \
$(RTRAC)/ve_gas_rt.o \
$(RTRAC)/wetdeprt.o \
$(RTRAC)/wrrcprt.o \
$(SOAP)/soachem.o \
$(SOAP)/soap.o \
$(SOAP)/spfcn.o \
$(SOAP)/vbschem.o \
$(HG)/hgadsorb.o \
$(HG)/hggaschem.o \
$(HG)/hgaqschem.o \
$(CAMX_MPI)/after_massum.o \
$(CAMX_MPI)/before_massum.o \
$(CAMX_MPI)/conc_update.o \
$(CAMX_MPI)/edge_pass.o \
$(CAMX_MPI)/find_edge.o \
$(CAMX_MPI)/init_lslice.o \
$(CAMX_MPI)/isbound.o \
$(CAMX_MPI)/master_bndconc.o \
$(CAMX_MPI)/master_send_1species_data.o \
$(CAMX_MPI)/master_send_gridded_data.o \
$(CAMX_MPI)/master_update.o \
$(CAMX_MPI)/newgrid.o \
$(CAMX_MPI)/node_get_feed.o \
$(CAMX_MPI)/node_get_lbc.o \
$(CAMX_MPI)/node_get_nbc.o \
$(CAMX_MPI)/node_recv_1species_data.o \
$(CAMX_MPI)/node_recv_gridded_data.o \
$(CAMX_MPI)/nodes_o3col.o \
$(CAMX_MPI)/nodes_alloc.o \
$(CAMX_MPI)/nodes_bndconc.o \
$(CAMX_MPI)/nodes_emiss.o \
$(CAMX_MPI)/node_send_feed.o \
$(CAMX_MPI)/node_send_lbc.o \
$(CAMX_MPI)/node_send_nbc.o \
$(CAMX_MPI)/nodes_init.o \
$(CAMX_MPI)/nodes_ipr.update.o \
$(CAMX_MPI)/nodes_irr.update.o \
$(CAMX_MPI)/nodes_met.o \
$(CAMX_MPI)/nodes_met_pig.o \
$(CAMX_MPI)/nodes_pass.o \
$(CAMX_MPI)/nodes_pass_sapnts.o \
$(CAMX_MPI)/nodes_pig_pass.o \
$(CAMX_MPI)/nodes_send_pig_sample.o \
$(CAMX_MPI)/nodes_send_rt_sample.o \
$(CAMX_MPI)/nodes_send_rtrcp_back.o \
$(CAMX_MPI)/nodes_send_walls_back.o \
$(CAMX_MPI)/nodes_topc.o \
$(CAMX_MPI)/nodes_tstep.o \
$(CAMX_MPI)/nodemass.o \
$(CAMX_MPI)/broadcast_grid_dimens.o \
$(CAMX_MPI)/broadcast_procid.o \
$(CAMX_MPI)/domain_decomp.o \
$(CAMX_MPI)/iniptr_node.o \
$(CAMX_MPI)/init_fields_node.o \
$(CAMX_MPI)/master_pass_xmass.o \
$(CAMX_MPI)/master_recv_1species_data.o \
$(CAMX_MPI)/myoffset.o \
$(CAMX_MPI)/node_send_1species_data.o \
$(CAMX_MPI)/nodes_send_lslice_back.o \
$(CAMX_MPI)/nodes_send_pig_back.o \
$(CAMX_MPI)/nodes_send_pig_misc.o \
$(CAMX_MPI)/nodes_send_pig_misc_real8.o \
$(CAMX_MPI)/nodes_send_pig_sum.o \
$(CAMX_MPI)/nodes_send_puffmass_back.o \
$(CAMX_MPI)/nodes_pass1.o \
$(CAMX_MPI)/mpi_feedback.o \
$(CAMX_MPI)/par_decomp.o \
$(CAMX_MPI)/par_decomp_bounds.o


model:  $(OBJCTS)
ifeq ($(MPI),false)
	$(CC) -c -o $(DUM_MPI)/mpi_dummy_c.o $(DUM_MPI)/mpi_dummy_c.c
	$(FC) -c $(FLGS) -o $(DUM_MPI)/mpi_dummy.o $(DUM_MPI)/mpi_dummy.f
else
	cd MPI/util;  make MPI_INST=$(MPI_INST) CC=$(CC); cd ../..
endif
ifeq ($(NCF),false)
	$(FC) -c $(FLGS) -o $(DUM_NCF)/ncf_dummy.o $(DUM_NCF)/ncf_dummy.f
endif
	$(FC) -o $(TARGT) $(FLGS) $(FLGS_OMP) $(FLGS_IEEE) $(MOD_OBJCTS) $(OBJCTS) $(LIBS) $(MPI_DUM_F) $(MPI_DUM_C) $(NCF_DUM_F)

modules:  $(MOD_OBJCTS)

.F90.o  :
	$(FC) -c -o $@ $(FLGS) $(CHUNK) $(FLGS_OMP) $(FLGS_IEEE) $(MODULES) $<
.f90.o  :
	$(FC) -c -o $@ $(FLGS) $(FLGS_OMP) $(FLGS_IEEE) $(MODULES) $<
.f.o    :
	$(FC) -c -o $@ $(FLGS) $(FLGS_OMP) $(FLGS_IEEE) $(MODULES) $<
.c.o    :
	$(CC) -c -o $@ $(INCLUDES) $<

