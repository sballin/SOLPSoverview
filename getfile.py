#!/bin/tcsh -f
# Extract the specified file from a case symlinked in ./history and place in current directory, possibly overwriting.
# Usage: in run directory, do getfile <filename e.g. b2mn.dat> <archive number e.g. 5>
tar -zxvf history/$2/20*.tar.gz ./$1
