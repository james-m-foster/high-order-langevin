Numerical simulation of the underdamped Langevin diffusion

This code estimates the discretization error of a new high order method for the underdamped Langevin equation:

dQ_t = P_t dt,
dP_t = -grad(U)(Q_t) dt - nu P_t dt + sqrt(2*nu / beta) dW_t.


The numerical method and example are outlined in the document langevin_presentation.pdf.

The source file langevin.cpp only requires headers from the C++ standard library.

The text file langevin_simulation.txt displays the output of the code.


License

The content of langevin_presentation.pdf is licensed under the Creative Commons Attribution 4.0 International license, and langevin.cpp is licensed under the MIT license.