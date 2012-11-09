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
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
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
	${OBJECTDIR}/main_haplo.o \
	${OBJECTDIR}/file_handling.o \
	${OBJECTDIR}/error.o \
	${OBJECTDIR}/ld.o


# C Compiler Flags
CFLAGS=-lcurl -Wl,-Bsymbolic-functions -lcprops -fopenmp -lm -lxml2 -std=gnu99 -largtable2 -lconfig

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=../../libs/common-libs/commons/file_utils.o ../../libs/common-libs/commons/http_utils.o ../../libs/common-libs/commons/result.o ../../libs/common-libs/commons/string_utils.o ../../libs/common-libs/commons/system_utils.o ../../libs/common-libs/containers/array_list.o ../../libs/common-libs/containers/list.o ../../libs/common-libs/containers/log.o ../../libs/bioinfo-libs/bioformats/vcf/vcf_file.o ../../libs/bioinfo-libs/bioformats/vcf/vcf_file_structure.o ../../libs/bioinfo-libs/bioformats/vcf/vcf_filters.o ../../libs/bioinfo-libs/bioformats/vcf/vcf_ragel.o ../../libs/bioinfo-libs/bioformats/vcf/vcf_reader.o ../../libs/bioinfo-libs/bioformats/vcf/vcf_stats.o ../../libs/bioinfo-libs/bioformats/vcf/vcf_util.o ../../libs/bioinfo-libs/bioformats/vcf/vcf_write.o ../../libs/bioinfo-libs/bioformats/features/region/region.o ../../libs/bioinfo-libs/bioformats/features/region/region_table.o ../../libs/bioinfo-libs/bioformats/features/region/region_table_utils.o ../../libs/bioinfo-libs/bioformats/gff/gff_batch.o ../../libs/bioinfo-libs/bioformats/gff/gff_file.o ../../libs/bioinfo-libs/bioformats/gff/gff_read.o ../../libs/bioinfo-libs/bioformats/gff/gff_reader.o ../../libs/bioinfo-libs/bioformats/gff/gff_write.o -lcprops -lconfig -largtable2 ../../global_options.o ../../hpg_variant_utils.o

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/common-libs/commons/file_utils.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/common-libs/commons/http_utils.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/common-libs/commons/result.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/common-libs/commons/string_utils.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/common-libs/commons/system_utils.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/common-libs/containers/array_list.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/common-libs/containers/list.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/common-libs/containers/log.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/vcf/vcf_file.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/vcf/vcf_file_structure.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/vcf/vcf_filters.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/vcf/vcf_ragel.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/vcf/vcf_reader.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/vcf/vcf_stats.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/vcf/vcf_util.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/vcf/vcf_write.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/features/region/region.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/features/region/region_table.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/features/region/region_table_utils.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/gff/gff_batch.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/gff/gff_file.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/gff/gff_read.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/gff/gff_reader.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../libs/bioinfo-libs/bioformats/gff/gff_write.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../global_options.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ../../hpg_variant_utils.o

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.c} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/main_haplo.o: main_haplo.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.c) -g -I../../libs/bioinfo-libs -I../../libs/common-libs -I.. -I/usr/local/include -I/usr/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main_haplo.o main_haplo.c

${OBJECTDIR}/file_handling.o: file_handling.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.c) -g -I../../libs/bioinfo-libs -I../../libs/common-libs -I.. -I/usr/local/include -I/usr/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/file_handling.o file_handling.c

${OBJECTDIR}/error.o: error.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.c) -g -I../../libs/bioinfo-libs -I../../libs/common-libs -I.. -I/usr/local/include -I/usr/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/error.o error.c

${OBJECTDIR}/ld.o: ld.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.c) -g -I../../libs/bioinfo-libs -I../../libs/common-libs -I.. -I/usr/local/include -I/usr/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ld.o ld.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/haplo

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
