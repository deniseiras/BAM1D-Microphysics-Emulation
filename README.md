# Machine Learning to Emulate Morrisson Cloud Microphysics implemented in BAM1D

This work aims to model an Artificial Neural Network, using Keras framework, for emulating the Morrisson Cloud Microphysics implemented in BAM1D (Brazilian Atmosphere Model - 1d column).

Some efforts where made in advance of emulating with lower errors using Neural Networks, as MLP, CNN1D, CNN2D, LSTM 

A version of BAM1D implements the emulation written in Fortran code, through the [Fortran Keras Bridge(FKB)](https://github.com/scientific-computing/FKB), which aims to be fast, for production, but work need to be done to accelerate the code using OpenMP.

There is another version, for research, implemented in Python, which communicates to BAM1D Fortran code through txt files. This version is more flexible, because the Python code receives the inputs from Fortran, calculates the results using any of the models (MLP, CNN, LSTM) and then generates an output file, which is read from Fortran file.

TODO:
- This text would be improved to resume better the applications
- Upload the BAM1D versions 
- Upload the paper after some corrections

<!-- This notebook has been saved in gdrive and copied here  -->
