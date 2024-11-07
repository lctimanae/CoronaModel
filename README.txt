If the code is used for a scientific work, please reference the following publication:

[1] Luis Carlos Timaná, Jaimis S. L. Colqui, Carlos Suárez, José Pissolato Filho. Implementation of a Voltage-Dependent Transmission Line Model with Corona Effect Consideration. Electric Power Systems Research (submitted article).

The transmission line models described in [1] are validated using high-voltage impulse test data presented in:

[2] C. F. Wagner, I. W. Gross, and B. L. Lloyd, "High-Voltage Impulse Tests on Transmission Lines [includes discussion]," in Transactions of the American Institute of Electrical Engineers. Part III: Power Apparatus and Systems, vol. 73, no. 2, pp. 196-210, April 1954.

The model in Section 2 of [1] is an improved version of the model presented in:

[3] T. M. Pereira and M. C. Tavares, "Development of a Voltage-Dependent Line Model to Represent the Corona Effect in Electromagnetic Transient Program," in IEEE Transactions on Power Delivery, vol. 36, no. 2, pp. 731-739, April 2021.

-----------------------------------------------------------------------------------------

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









