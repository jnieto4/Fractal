# Fractal
This program uses MPI parallel programming to create a fractal image

Steps to run (using debian based systems):

sudo apt update
sudo apt install openmpi-bin libopenmpi-dev g++

mpic++ fractal_MPI.cpp -o fractal
mpirun -np 4 ./fractal x


x is any value between 0-1024, this determines the width/height in pixles the result will be
