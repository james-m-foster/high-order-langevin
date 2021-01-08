Numerical simulation of the underdamped Langevin diffusion

dx_t = x_t dt,
dv_t = - gamma v_t dt - u grad(f)(x_t) dt + sqrt(2\gamma u) dW_t.

The R scripts in the logistic_regresssion folder reproduce the results of the paper "The shifted ODE method for underdamped Langevin MCMC" by James Foster, Terry Lyons and Harald Oberhauser.

The scripts SOFA_method.R, Strang_splitting.R, Randomized_midpoint_method.R, OBABO_scheme.R and Leftpoint_method.R generate and compare sample paths of logistic regression weights (for each method). The data is loaded by GermanCredit.R using the package unbiasedmcmc.


The code in the thesis_experiment folder reproduces the double-well example in Section 5.3 of my DPhil thesis.

The source file langevin.cpp only requires headers from the C++ standard library.

The text file langevin_simulation.txt displays the output of langevin.cpp.

The numerical methods and example for the thesis are outlined in presentation/langevin_presentation.pdf.


References

The unbiasedmcmc package is based on the paper "Unbiased Markov chain Monte Carlo with couplings", by Pierre E. Jacob, John O'Leary and Yves F. Atchade.


License

The content of langevin_presentation.pdf is licensed under the Creative Commons Attribution 4.0 International license and OFA_method.R, Strang_splitting.R, Randomized_midpoint_method.R, OBABO_scheme.R, Leftpoint_method.R, GermanCredit.R and langevin.cpp are licensed under the MIT license.