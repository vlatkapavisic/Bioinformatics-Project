Hirschberg.exe contains code for local alignment with sampling a pair of sequences
using Smith-Waterman algorithm, followed by Hirschberg algorithm using Needleman-Wunsch algorithm.

There are a few options for input:
1. two file paths, containing one sequence each, either as plain text, or in FASTA format
2. same as 1., with third file path on the end, which represents the output file path

optional parameters:
"-c" causes output alignment to console, along with score and matching
"-t" causes measuring of time, and outputs measured miliseconds in console, regardless of -c 
"-m" causes measuring of memory, and outputs measured bytes in console, regardless of -c

Measuring memory will affect time
Alignment output will severely affect time for long strings

Input file paths can be absolute or relative,
Output file path can be absolute or relative,
both regardless of file extension, but extension must 
be entered as a part of filename in arguments while starting the program

checking mono version: mono --version
installing mono: sudo apt-get install monodevelop 
(installed on Ubuntu 12.04, and bioLinux from http://nebc.nerc.ac.uk/tools/bio-linux/)
compiling: dmcs filename.cs -unsafe
running: 
mono filename.cs inputFilePath1 inputFilePath2 
or
mono filename.cs inputFilePath1 inputFilePath2 outputFilePath