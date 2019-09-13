#include <fstream>
#include <iomanip>
#include <random>

using namespace std;

// Numerical simulation of the underdamped Langevin diffusion on [0, T]:
// dP_t = Q_t dt
// dQ_t = -grad(U)(Q_t) dt - nu P_t dt + sqrt(2*nu / beta) dW_t

class ULangevinMethods {
    public:
        // Input parameters
        double nu, beta, T;
        int no_of_steps;

        ULangevinMethods(double, double, double, int );
        double grad_potential(double);
        double ode_vect_field(double, double, double, double, double);
        pair<double, double> linear_ode_step(pair<double, double>, double, bool);
        pair<double, double> zigzag_ode_step(pair<double, double>, double, double, double);

    private:
        // High order piecewise linear path coefficients
        const double a = 0.5 - sqrt(0.05);
        const double one_minus_a = 0.5 + sqrt(0.05);
        const double one_minus_two_a = sqrt(0.2);
        const double half_of_one_minus_two_a = sqrt(0.05);
        const double one_over_one_minus_a = 2.5 - sqrt(1.25);
        const double six_over_one_minus_a = 15.0 - sqrt(45.0);

        // Runge-Kutta-Nystrom coefficients
        const double a21 = 0.164217030;
        const double a31 = 0.139710559;
        const double a32 = -0.103005664;
        const double b1 = 0.260311692;
        const double b2 = 0.4039053382;
        const double b3 = -0.164217030;
        const double c2 = 0.630847693;
        const double c3 = 0.536704894;
        const double d1 = 0.260311692;
        const double d2 = 1.094142798;
        const double d3 = -0.354454490;

        // Precomputed values that depend on the input parameters
        double sigma, quarter_nu_squared, half_nu;
        double step_size;

        // Precomputed values that depend on step_size
        double small_step_size;
        double small_step_size_a21;
        double small_step_size_a31;
        double small_step_size_a32;
        double small_step_size_c2;
        double small_step_size_c3;
        double exp_minus_small;
        double exp_minus_small_c2;
        double exp_small_c2;
        double exp_minus_small_c3;
        double exp_small_c3;

        double big_step_size;
        double big_step_size_a21;
        double big_step_size_a31;
        double big_step_size_a32;
        double big_step_size_c2;
        double big_step_size_c3;
        double exp_minus_big;
        double exp_minus_big_c2;
        double exp_big_c2;
        double exp_minus_big_c3;
        double exp_big_c3;
};

// Constructor will compute the above private variables
ULangevinMethods::ULangevinMethods(double input_nu, double input_beta,
                                   double input_T, int input_no_of_steps){
    nu = input_nu;
    beta = input_beta;
    T = input_T;
    no_of_steps = input_no_of_steps;

    sigma = sqrt(2.0*input_nu/input_beta);
    quarter_nu_squared = 0.25*pow(input_nu, 2);
    half_nu = 0.5*input_nu;

    step_size =  input_T/(double)input_no_of_steps;

    // The "change of variable" parameters
    small_step_size = a*step_size;
    small_step_size_a21 = small_step_size*a21;
    small_step_size_a31 = small_step_size*a31;
    small_step_size_a32 = small_step_size*a32;
    small_step_size_c2 = small_step_size*c2;
    small_step_size_c3 = small_step_size*c3;
    exp_minus_small = exp(-half_nu*small_step_size);
    exp_minus_small_c2 = exp(-half_nu*small_step_size_c2);
    exp_small_c2 = exp(half_nu*small_step_size_c2);
    exp_minus_small_c3 = exp(-half_nu*small_step_size_c3);
    exp_small_c3 = exp(half_nu*small_step_size_c3);

    big_step_size = one_minus_two_a*step_size;
    big_step_size_a21 = big_step_size*a21;
    big_step_size_a31 = big_step_size*a31;
    big_step_size_a32 = big_step_size*a32;
    big_step_size_c2 = big_step_size*c2;
    big_step_size_c3 = big_step_size*c3;
    exp_minus_big = exp(-half_nu*big_step_size);
    exp_minus_big_c2 = exp(-half_nu*big_step_size_c2);
    exp_big_c2 = exp(half_nu*big_step_size_c2);
    exp_minus_big_c3 = exp(-half_nu*big_step_size_c3);
    exp_big_c3 = exp(half_nu*big_step_size_c3);
}

// Gradient of a double-well potential U(q) = (q^2 -1)^2
double ULangevinMethods::grad_potential(double x){
    return 4.0*x*(pow(x,2) - 1.0);
}

// Vector field of the second order ODE (i.e. after the change of variables)
double ULangevinMethods::ode_vect_field(double y,
                                        double time_increment, double brownian_increment,
                                        double exp_minus_nu, double exp_nu){

    return (quarter_nu_squared*y - exp_nu*grad_potential(exp_minus_nu*y)) \
                 *time_increment + exp_nu*sigma*brownian_increment;
}

// Method for propagating the numerical solution along a linear path
pair<double, double> ULangevinMethods::linear_ode_step(pair<double, double> qp,
                                                       double path_increment,
                                                       bool small_step){

    // Use change of variables to get initial values for the ODE: y" = F(y)
    double y = qp.first;
    double y_prime = half_nu*qp.first + qp.second;

    // Apply a step of the third order symplectic RKN method
    // Note this is only performed for two sizes of interval
    if (small_step){
        double k1 = ode_vect_field(y, small_step_size, path_increment, 1.0, 1.0);

        double captial_y_2 = y + small_step_size_c2*y_prime + small_step_size_a21*k1;

        double k2 = ode_vect_field(captial_y_2,
                                   small_step_size, path_increment,
                                   exp_minus_small_c2, exp_small_c2);

        double captial_y_3 = y + small_step_size_c3*y_prime \
                               + small_step_size_a31*k1 + small_step_size_a32*k2;

        double k3 = ode_vect_field(captial_y_3,
                                   small_step_size, path_increment,
                                   exp_minus_small_c3, exp_small_c3);

        y = y + small_step_size*(y_prime + b1*k1 + b2*k2 + b3*k3);
        y_prime = y_prime + d1*k1 + d2*k2 + d3*k3;

        // Use change of variables to get back (q, p)
        y = exp_minus_small*y;
        y_prime = exp_minus_small*y_prime - half_nu*y;
    }
    else
    {
        double k1 = ode_vect_field(y, big_step_size, path_increment, 1.0, 1.0);

        double captial_y_2 = y + big_step_size_c2*y_prime + big_step_size_a21*k1;

        double k2 = ode_vect_field(captial_y_2,
                                   big_step_size, path_increment,
                                   exp_minus_big_c2, exp_big_c2);

        double captial_y_3 = y + big_step_size_c3*y_prime \
                               + big_step_size_a31*k1 + big_step_size_a32*k2;

        double k3 = ode_vect_field(captial_y_3, big_step_size, path_increment,
                                   exp_minus_big_c3, exp_big_c3);

        y = y + big_step_size*(y_prime + b1*k1 + b2*k2 + b3*k3);
        y_prime = y_prime + d1*k1 + d2*k2 + d3*k3;

        // Use change of variables to get back (q, p)
        y = exp_minus_big*y;
        y_prime = exp_minus_big*y_prime - half_nu*y;
    }

    return make_pair(y, y_prime);
};

// Method for propagating the numerical solution along each part of \hat{W}
pair<double, double> ULangevinMethods::zigzag_ode_step(pair<double, double>  qp,
                                                       double brownian_increment,
                                                       double brownian_area,
                                                       double brownian_skew_area){

    // Compute the connecting points for the discretized Brownian path \hat{W}
    double half_of_cplusb = 0.5*brownian_increment \
                              + one_over_one_minus_a*brownian_area;

    double half_of_cminusb = half_of_one_minus_two_a*brownian_increment \
                              - six_over_one_minus_a*brownian_skew_area;

    double b = half_of_cplusb - half_of_cminusb;
    double c = half_of_cplusb + half_of_cminusb;

    pair<double, double> y = qp;

    // Propagate the numerical solution along the piecewise linear path
    y = linear_ode_step(y, b, true);
    y = linear_ode_step(y, c - b, false);
    y = linear_ode_step(y, brownian_increment - c, true);

    return y;
};

int main()
{
    // Input parameters
    const double nu = 1.0;
    const double beta = 3.0;
    const pair<double, double> y0 = make_pair(0.0, 0.0);
    const double T = 10.0;
    const int no_of_steps = 1000;

    // Number of steps used by the "fine" approximation
    // during each step of the "crude" numerical method
    const int no_of_fine_steps = 50;

    // Number of paths used for the Monte Carlo estimators
    const int no_of_paths = 100000;

    // Variances for generating the normal random variables
    const double third = 1.0/3.0;
    const double twelve = 1.0/12.0;
    const double seven_twentieth = 1.0/720.0;

    // Step size parameters
    const double step_size =  T/(double)no_of_steps;
    const double one_over_step_size = 1.0/step_size;
    const double one_over_step_size_squared = pow(one_over_step_size, 2);
    const double fine_step_size = T/(double)(no_of_steps*no_of_fine_steps);
    const double half_fine_step_size = 0.5*fine_step_size;
    const double fine_step_size_squared =  pow(fine_step_size, 2);

    // We will be comparing the methods on two different time scales
    ULangevinMethods course_method(nu, beta, T, no_of_steps);
    ULangevinMethods fine_method(nu, beta, T, no_of_steps*no_of_fine_steps);

    // Normal distributions for generating the various increments and time areas
    std::random_device rd;
    std::default_random_engine generator(rd());

    std::normal_distribution<double> increment_distribution(0.0,sqrt(step_size));
    std::normal_distribution<double> area_distribution(0.0, sqrt(twelve*step_size));
    std::normal_distribution<double> skew_area_distribution(0.0, sqrt(seven_twentieth \
                                                                       *step_size));

    std::normal_distribution<double> fine_increment_distribution(0.0, sqrt(fine_step_size));
    std::normal_distribution<double> fine_area_distribution(0.0, sqrt(twelve*fine_step_size));
    std::normal_distribution<double> fine_skew_area_distribution(0.0, sqrt(seven_twentieth \
                                                                            *fine_step_size));

    // Numerical solutions computed with course and fine step sizes
    pair<double, double> y_ODE = y0;
    pair<double, double> y_fine = y0;

    // Information about the Brownian motion (increments and areas)
    // These objects are described in langevin_presentation.pdf as:
    // brownian_increment is W_{s,t}
    // brownian_area      is H_{s,t}
    // brownian_skew_area is K_{s,t}
    double brownian_increment = 0.0;
    double brownian_area = 0.0;
    double brownian_skew_area = 0.0;
    double fine_brownian_increment = 0.0;
    double fine_brownian_area = 0.0;
    double fine_brownian_skew_area = 0.0;

    // Time variable used in each "fine" simulation
    double fine_time = 0.0;

    // Strong and weak error estimators for y_ODE at time T
    double end_point_error = 0.0;
    double sec_moment_error = 0.0;

    double samplepath[no_of_steps + 1] = {0.0};

    for (int i=0; i<no_of_paths; ++i) {
        for (int j=1; j<=no_of_steps; ++j) {

            brownian_increment = 0.0;
            brownian_area = 0.0;
            brownian_skew_area = 0.0;
            fine_time = 0.0;

            for (int k=0; k < no_of_fine_steps; ++k){
                // Generate information about the Brownian path over the "fine" increment
                fine_brownian_increment = fine_increment_distribution(generator);
                fine_brownian_area = fine_area_distribution(generator);
                fine_brownian_skew_area = fine_skew_area_distribution(generator);

                // Propagate the numerical solution over the fine increment
                y_fine = fine_method.zigzag_ode_step(y_fine, fine_brownian_increment,
                                                     fine_brownian_area, fine_brownian_skew_area);

                // Update the information about the Brownian path over the
                // course increment using the recently generated variables.
                // The below procedure can be derived using some elementary
                // properties of integration (additivity and linearity)
                brownian_skew_area = brownian_skew_area \
                                   + fine_step_size_squared \
                                      * (third*fine_brownian_increment \
                                         + 0.5*fine_brownian_area - fine_brownian_skew_area) \
                                   + fine_step_size \
                                      * (half_fine_step_size*brownian_increment + fine_time \
                                         *(brownian_increment \
                                            + 0.5*fine_brownian_increment + fine_brownian_area));

                brownian_area = brownian_area + fine_step_size \
                                 * (brownian_increment \
                                     + 0.5*fine_brownian_increment + fine_brownian_area);

                brownian_increment = brownian_increment + fine_brownian_increment;

                fine_time = fine_time + fine_step_size;
            }

            // Compute the time areas for the Brownian path over the course increment
            brownian_area = brownian_area*one_over_step_size - 0.5*brownian_increment;

            brownian_skew_area = third*brownian_increment + 0.5*brownian_area \
                                  - brownian_skew_area*one_over_step_size_squared;

            // Propagate the numerical solution over the course increment
            y_ODE = course_method.zigzag_ode_step(y_ODE, brownian_increment,
                                                  brownian_area, brownian_skew_area);

            // Store the sample path if we have reached the final iteration
            if (i == no_of_paths - 1){
                samplepath[j] = y_ODE.first;
            }
        }

        // Update the error estimators
        end_point_error = end_point_error + pow(y_ODE.first - y_fine.first, 2);
        sec_moment_error = sec_moment_error + pow(y_ODE.first, 2) - pow(y_fine.first, 2);

        // Reset the numerical solutions
        y_ODE = y0;
        y_fine = y0;
    }

    // Compute the strong and weak error estimators
    end_point_error = sqrt(end_point_error/(double)no_of_paths);
    sec_moment_error = abs(sec_moment_error/(double)no_of_paths);

    // Display the results in a text file
    ofstream myfile;
    myfile.open ("langevin_simulation.txt");

    myfile << std::fixed << std::setprecision(2) \
           << "L2 error at time T = " << T << ": \t\t\t"
           << std::setprecision(15) << end_point_error << "\n";

    myfile << std::fixed << std::setprecision(2) \
           << "Error in second moments at time T = " << T << ": \t"
           << std::setprecision(15) << sec_moment_error << "\n";

    myfile << std::fixed << std::setprecision(15) \
           << "Number of steps: " << "\t\t\t\t" << no_of_steps << "\n\n";

    myfile << std::fixed << std::setprecision(15) \
           << "Example sample path" << "\n\n" ;

    myfile << std::fixed << std::setprecision(15) \
            << "t" << "\t\t\t" << "Q_t" << "\n";

    for (int j=0; j<=no_of_steps; ++j) {
        myfile << std::fixed << std::setprecision(15) \
               << j*step_size << "\t" << samplepath[j] << "\n";
    }

    myfile.close();

    return 0;
}