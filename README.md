# missile_estimation

## Abstract
The goal of this project was to intercept a target with a missile using line-of sight mea-
surements. Two models were used to develop the state dynamics in this project: the Gauss-
Markov model and Random Telegraph model. The Continuous-time Kalman Filter was used
to determine minimum variance estimates of the lateral position, velocity, and target accel-
eration for both models. These estimates were compared to their corresponding true values
over a span of 10 seconds. A Monte Carlo simulation was run for 10,000 realizations in order
to confirm the Kalman Filter algorithm was functional and the models used were approx-
imately correct. A comparison of the root mean square error of the different states from
the simulation is made with the corresponding filter values. The Kalman Filter performed
similarly for simulations of both the Gauss-Markov model and Random Telegraph model.

## Code Descriptions

**one_real_gm** runs one realization of the continuous time kalman filter using the gauss-markov model. It plots the true states vs. the estimated states

**one_real_tele** runs one realization of the continuous time kalman filter using the random telegraph model. It plots the true states vs. the estimated states

**monte_carlo_gm** runs a monte carlo simulation of the continuous time kalman filter using the gauss-markov model. It plots the simulated root mean square error of the states versus those calculated by the filter.

**monte_carlo_tele** runs a monte carlo simulation of the continuous time kalman filter using the telegraph model. It plots the simulated root mean square error of the states versus those calculated by the filter.

**ct_kalman_filter** is a helper function used in monte_carlo_gm to automate calculations

**tele_kalman_filt** is a helper function used in monte_carlo_tele to automate calculations
