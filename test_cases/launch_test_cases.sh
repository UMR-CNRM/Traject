#!/user/bin/bash
#Full test cases for traject
#This file describes the test cases and processes them
#The logs and output files are generated in a dedicated directory (./tests_cases/test{N}/)

#PYTHONPATH=../src:$PYTHONPATH

########################
#Test case 1 : mslp algorithm on extratropical cyclone Alex on an ARPEGE +48h forecast
#########################
if [ 0 == 0 ]
then
repout=./test1/
echo "##############" > $repout/test_case1.log
echo "Running test case 1 " >> $repout/test_case1.log
echo "##############" >> $repout/test_case1.log
echo "Running test case 1 ... output files in " $repout
python3 test_case1.py >> $repout/test_case1.log
echo "Test case 1 OK "
echo "Test case 1 OK " >> $repout/test_case1.log
fi
#########################
#Test case 2 : VDG algorithm on extratropical cyclone Alex on two PEARP members at +48h forecast
#(try initial time : 20201001-06UTC : one member starts the track from this instant, one member starts later)
#otherwise 20201001-12UTC ...
#########################
if [ 0 == 0 ]
then
repout=./test2/
echo "##############" > $repout/test_case2.log
echo "Running test case 2 " >> $repout/test_case2.log
echo "##############" >> $repout/test_case2.log
echo "Running test case 2 ... output files in " $repout
python3 test_case2.py >> $repout/test_case2.log
echo "Test case 2 OK "
echo "Test case 2 OK " >> $repout/test_case2.log
fi
#########################
#Test case 3 : Ayrault algorithm in ARPEGE analyses: we should find the ALEX track (loop aloft Brittany)
#########################
if [ 0 == 0 ]
then
repout=./test3/
echo "##############" > $repout/test_case3.log
echo "Running test case 3 " >> $repout/test_case3.log
echo "##############" >> $repout/test_case3.log
echo "Running test case 3 ... output files in " $repout
python3 test_case3.py >> $repout/test_case3.log
echo "Test case 3 OK "
echo "Test case 3 OK " >> $repout/test_case3.log
fi
#########################
#Test case 4 : Ayrault algorithm January 2020 (Mediterranean cyclones) on ERA5 data
#########################
if [ 0 == 0 ]
then
repout=./test4/
echo "##############" > $repout/test_case4.log
echo "Running test case 4 " >> $repout/test_case4.log
echo "##############" >> $repout/test_case4.log
echo "Running test case 4 ... output files in " $repout
python3 test_case4.py >> $repout/test_case4.log
echo "Test case 4 OK "
echo "Test case 4 OK " >> $repout/test_case4.log
fi

