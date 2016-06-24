#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=nvcc
CCC=nvcc
CXX=nvcc
FC=gfortran
AS=as

# Macros
CND_PLATFORM=CUDA-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/GraphHNMain.o \
	${OBJECTDIR}/GraphHullNumber.o \
	${OBJECTDIR}/UndirectedSparseGraph.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/graph-hull-number-parallel

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/graph-hull-number-parallel: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/graph-hull-number-parallel ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/GraphHNMain.o: GraphHNMain.cu 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/GraphHNMain.o GraphHNMain.cu

${OBJECTDIR}/GraphHullNumber.o: GraphHullNumber.cu 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.c) -g -o ${OBJECTDIR}/GraphHullNumber.o GraphHullNumber.cu

${OBJECTDIR}/UndirectedSparseGraph.o: UndirectedSparseGraph.cpp 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.cc) -g -o ${OBJECTDIR}/UndirectedSparseGraph.o UndirectedSparseGraph.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/graph-hull-number-parallel

# Subprojects
.clean-subprojects:
