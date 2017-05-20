#!/bin/sh
#
##set -x
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
# Get also the absolute path of the object directory
#
objdir=$(
cd -P -- "$(pwd)" &&
pwd -P
)
#
destdir=$objdir
#
# Replicate the hierarchy of makefiles for the library only
#
(cd $srcdir;
   for i in $(find .  \
		 -name \[mM\]akefile | grep -v \\./Makefile); do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp $relpath/*akefile ${destdir}/$relpath
   done
)
#
# Copy full Examples material
#
(cd $srcdir; cd .. ; cp -rp examples ${destdir} )
#
# Copy other needed top-level files
#
cp -p $srcdir/psml.mk.in ${destdir}
#
# Set the appropiate variables in the build makefile
#
sed "s#VPATH=\.#VPATH=${srcdir}#g" ${srcdir}/makefile | \
sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${destdir}/makefile

#
echo " *** Compilation setup done. "

## for i in $(find . \( -path Tutorial -o -path Examples \) -prune -o \
