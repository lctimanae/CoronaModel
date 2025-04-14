The main files and their functions are described below:

Test.m: Simulates high-voltage impulse tests on transmission lines as described in [2], using the models presented in [1].

exper_simul_comparison_ini.m: Simulates a transmission line with its ends connected to a matching resistor and a voltage source representing overvoltage caused by a lightning surge. Uses the VDLM model described in [1].

exper_simul_comparison_2.m: Simulates a transmission line with its terminals connected to a matching resistor and a voltage source representing overvoltage caused by lightning. Uses the AVDLM model presented in [1].

BergeronModel_ini.m: VDLM model for simulating transmission lines considering the corona effect [1].

BergeronModel.m: AVDLM model for simulating transmission lines considering the corona effect [1].

error_calculation.m: Calculates the mean percentage error between simulated and experimental values.

Results/SimulationResults.m: Shows a comparison of simulation results from the VDLM and AVDLM models presented in [1] with experimental data from [2].

Results/ResultsVDLM/VDLM_Comparison.m: Compares experimental results from [2] with simulations from the VDLM and AVDLM models presented in [1] and [3].

Results/ResultsVDLM/error_calculation.m: Calculates the mean percentage error between simulated and experimental values.

Optimization/OptimizationGA.m: Adjusts the line model parameters. In [1] the optimized parameters of the AVDLM model are also used for the VDLM model.









