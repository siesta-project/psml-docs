#!/bin/sh
#
#  Run 'show_psml' on a file in 'plot' mode to generate
#  a log file and data files in a 'test_output' directory.
#  Its contents can be compared with those of 'reference_output'
#  to check the installation.
#
#  This script must be run in the (Building Directory)/examples directory,
#  after compilation of the example programs (see the installation
#  instructions):
#
#            sh run_test.sh
#
echo
echo "Running installation test..."
echo
#
if [ -d test_output ]
then
    echo "Directory 'test_output' must be removed"
    exit
fi

mkdir -p test_output
cd test_output
../show_psml -p ../Ba.sc-ionic-siesta-vnl.psml > Ba.sc-ionic-siesta-vnl.show
cd ..
#
echo "Test output in 'test_output'"
echo "Compare to contents of 'reference_output'"
echo "You might want to use 'diff -rq test_output reference_output'"

