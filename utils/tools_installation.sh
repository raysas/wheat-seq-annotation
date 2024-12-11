#!/bin/bash

# -- creates a directory if it doesn't exist to put it in tehre and keep it clean
if [ ! -d tools/ ]; 
    then mkdir tools; 
fi 


# -------------------- Genome Browser --------------------
# <> Artemis Genome Browser (https://www.sanger.ac.uk/science/tools/artemis)
#  (can also install it through conda)

wget https://github.com/sanger-pathogens/Artemis/releases/download/v18.2.0/artemis-unix-release-18.2.0.tar.gz #downloads the tar file
if [ ! -d tools/ ]; then mkdir tools; fi #creates a directory if it doesn't exist to put it in tehre and keep it clean

mv artemis-unix-release-18.2.0.tar.gz tools/
tar zxf tools/artemis-unix-release-18.2.0.tar.gz

# running it:
# ``` 
artemis/art 
#or art if installed through conda
# ```
# or try opening a fasta example from terminal:  
# ```
artemis/art utils/example/thaliana_example.fasta
# ```
