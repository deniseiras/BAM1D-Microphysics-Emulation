# Machine Learning to Emulate Morrisson Cloud Microphysics implemented in BAM1D

## Resume

The parameterized microphysics of a numerical weather or climate prediction model is one of the parts of the system that requires the greatest computational load. Furthermore, the understanding of physical laws is not fully represented in the parameterizations. In this work, some Deep Learning architectures and techniques were applied in order to emulate a microphysics parameterization of two-moment clouds in a Single Column Model (1D). The best solution found presented minor errors when comparing the prediction with the original microphysics.

## History

This work created Artificial Neural Networks models for use in emulation of the Morrisson Cloud Microphysics implemented in BAM1D.
BAM1D is a Single Column Model which implements the same Microphysics implemented in BAM3D. It was written in Fortran code.

Some efforts where made in a Python Jupyter Notebook, in advance of finding the best Neural Networks arhitectures, as MLP, CNN1D, CNN2D, LSTM to emulate the microphysics. 
The data used was generated from input and output of BAM1D microphysics routine, and then was used to train the Neural Networks. The best models where chosen to emulate the microphysics in BAM1D.

One of the BAM1D versions implements the emulation in Fortran code by using the [Fortran Keras Bridge(FKB)](https://github.com/scientific-computing/FKB) framework, that implements Keras MLP in Fortran. This version aims to be fast, for production level, but work need to be done to accelerate the code using OpenMP, for example.

There is another version of BAM1D, for research, which communicates to a Python code (bam1d_physics_emulation.py) through txt files. This version is more flexible, because Python code receives the inputs from Fortran, calculates the results using any of the models (MLP, CNN, LSTM) and then generates the emulations as an output file, which is read from BAM1D code in turn.

TODO:
- Upload the BAM1D versions 
- Upload the paper after some corrections

<!-- This notebook has been saved in gdrive and copied here  -->
