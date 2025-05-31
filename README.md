# ABTestFDR
Supplementary code for Design of Bayesian A/B Tests Controlling False Discovery Rates and Power manuscript 

There are five code files for this manuscript.

<<<<<<< HEAD
    01a_build_tree: creates a tree to characterize the admissible combinations of binary outcomes similar to Figure 1
    01b_lp: these functions can be used to obtain the multinomial models that define the data generation process
	01c_alg2: the functions to implement Algorithm 2 in various contexts
	02a_get_mixture: the code to get the 30 multinomial models that define our data generation process
	03_sim_sec5: the code to reproduce the numerical results in Section 5
=======
- 01a_build_tree: creates a tree to characterize the admissible combinations of binary outcomes similar to Figure 1
- 01b_lp: these functions can be used to obtain the multinomial models that define the data generation process
- 01c_alg2: the functions to implement Algorithm 2 in various contexts and construct bootstrap confidence intervals
- 02a_get_mixture: the code to get the 30 multinomial models that define our data generation process
- 03_sim_sec5: the code to reproduce the numerical results in Section 5
- 04_contour: the code to reproduce the contour plots in Section 6
>>>>>>> 1dafd98aa9c052868e8f335addd3b24183661ad7

There is also the "tree_opt.csv" file that characterizes the admissible outcomes for the Optimizely example and the "params_full.csv" file that details the multinomial parameters for all 30 models in our data generation process.
