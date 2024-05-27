#!/bin/sh

#
# This is not completely kosher. It will find the PSML file
# only if the "build" directory is a first-level subdirectory of the package directory.
# ToDo: use the proper auto-xxx variables.
# ToDo: Check that the output is really '83'
#
./getz ../../examples/83_Bi_r.psml
