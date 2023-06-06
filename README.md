# Multi-objective-COBYLA-optimization-for-supply-chain-mangment
it involves the implementation of pareto forntier for 3 objectives in the supply chain optimization 
the code represents the applicaiton of a multiobjective cobyla in supply chain . follow the guided steps in the note book for 
adding constraints and modeling the optimization problem.



## Background information about the stochastic sustainable supply chain system 
The closed-loop supply chain (CLSC) is the supply chain that includes various recovery plans for used products to be reused in the industry. Most of the previous stochastic CLSC studies considered the effect of uncertain parameter changes on the economic aspect only, while the other sustainability aspects were neglected. The purpose of this code is to develop a realistic mathematical model that represents and analyzes the impact of uncertainty in demand and recovery rate of products on the economic, environmental, and social sustainability aspects in the CLSC. The objective functions were optimized using the constrained optimization by linear approximation (COBYLA) algorithm along with preference-based Pareto optimal solution set algorithm at multiple computations to optimize various objectives simultaneously and efficiently. the CLSC stochastic model developed is represented as shown below in figure 1 and solved in figure 2






<img src="https://github.com/Omarelfarouk90/Multi-objective-COBYLA-optimization-for-supply-chain-mangment/assets/53394104/ee7b9ca0-a48e-4ad2-a895-969bd1a7efb8" width="600" height="200">


Figure 1: CLSC model describtion amont various tiers



<img src="https://github.com/Omarelfarouk90/Multi-objective-COBYLA-optimization-for-supply-chain-mangment/assets/53394104/8100e30b-4251-424b-ad1d-6a31e31dd641" width="500" height="500">


Figure 2: Flowchart  for solving the stocahstic CLSC model

The mathematical model developed for the stochastic CLSC and the Pareto COBYLA algorithm used are found in the journal with in depth details

Credits goes to Julian Black for developing the Pymoo code for verification of the Pareto COBYLA results with NSGA(Non- dominated sorting genetic algorithm)

### Link for the journal paper

https://www.tandfonline.com/doi/abs/10.1080/21681015.2021.1963338
