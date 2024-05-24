# Flow-based Offline Charging Scheduler (FOCS)

For large groups of electric vehicles (EVs), controllers may repeatedly solve the following offline problem:

Each EV charging session is a job j with arrival time, departure time, energy demand and maximum charging power. 
A feasible charging schedule has to:
1) schedule j only between its arrival and departure
2) schedule j to charge at least its energy demand between its arrival and departure
3) charge j at non-negative power (no V2G)
4) charge j at a rate no more than its maximum charging power
Furthermore, we want to make the aggregated power profile 'as flat as possible'. This may be the square of the l2-norm of the aggregated power (in the discretized setting). But it may also be any other separable function of the aggregated power that is convex, differentiable and increasing. 

FOCS [1] is an algorithm that determines an optimal solution for the problem sketched above. For more information, we refer you to the paper.
This repository holds a proof of concept implementation of FOCS [1]. 

## Code
The implementation of FOCS heavily relies on the python package networkx [2] for its network structures and maximum flow solvers. 
Before running it yourself, make sure to check for hard coded paths. We both read and write csv files and the code was build as a proof of concept, not for usability. 
Furthermore, make sure from the working directory of your code, you add 'data/input/', 'data/output/' and within the latter there should be subdirectories 'flowmethod', 'histograms', 'instances' and 'timesteps'. 

As input, the code requires a csv to read a pandas dataframe from. Each row represents a charging session, or more generally a job in the problem instance. The following headers should be present:
<table>
  <tr>
    <th>column header </th>
    <th>unit </th>
    <th>notes </th>
  </tr>

  <tr>
    <td>'total_energy_Wh' </td>
    <td>[Wh] </td>
    <td>Total Energy demand. Note that the FOCSinstance class automatically converts this to kWh. If numerical inputs are too large, that may result in errors.  </td>
  </tr>

  <tr>
    <td>'average_power_W' </td>
    <td>[W] </td>
    <td>Average power if energy demand served between arrival and departure time.  </td>
  </tr>

  <tr>
    <td>'maxPower' </td>
    <td>[W] </td>
    <td>Maximum power. FOCSinstance objects will compare the average and maximum input to decide on the actual maximum power.  </td>
  </tr>

  <tr>
    <td>'t0_[timestep]' </td>
    <td>[units of timeStep since start] </td>
    <td>Arrival time. Here, timeStep is the time granularity considered. For example, 900 would mean quarterly granularity, the column header would be 't0_900'. The arrival is at the beginning of the quarter. (Therefore including the quarter). </td>
  </tr>

  <tr>
    <td>'t1_[timestep]' </td>
    <td>[units of timeStep since start] </td>
    <td>Departure time. Here, timeStep is the time granularity considered. For example, 900 would mean quarterly granularity, the column header would be 't1_900'. The departure is at the beginning of the quarter. (Therefore excluding the quarter). </td>
  </tr>
</table>
		
In our initial experiments, we used data collected at an office location parking lot and processed it with the code shared here: https://github.com/lwinschermann/OfficeEVparkingLot. 
See dataGen.py for the code that called upon the other repository and the exact filtering settings.

The code is run from the main.py.
For comparison, the code includes a Gurobi [3] model that can solve the same problem. Note that it is misleadingly named lp.py. It is however not a linear program, but a quadratic one in this implementation!
secondary.py provides some code to combine individual experiments into batches. Furthermore, extensions AVR.py and OA.py contain online algorithms for the considered EV scheduling problem. Price.py provides an extension that incorporates CO2 objectives, both in itself and as a weighted combination with the flattening objective.

## Past experiments
### Runtime experiments
Experimental setup December 2023/January 2024 <br>
Run with Python 3.8.8 <br>
We conducted runtime experiments for the following parameters:

reps = 500<br>
instanceSizes = [n for n in range(1,20)] + [n for n in range(20,501,10)] + [n for n in range(600, 1001, 100)] <br>
timeSteps = [60,900,1800,3600]<br>
maxFlows = [shortest_augmenting_path, edmonds_karp, preflow_push, dinitz]<br>
write = True<br>
randomSample = True <br>

The various timeSteps and maxFlows were run in parallel to speed up copmletion time of all experiments. Under 'runtime_experiments/' we added the corresponding main files for those batches. To run them, you have to move them one repository up. <br>
The following CPU times were recorded (using time.process_time()) for experiments x_y, where x is the index of the maxFlow in the list above, and y is the timeSteps value. 

error margin 0.0000001<br>
0_3600 : 140754<br>
2_3600 : 147725<br>
1_3600 : 221068<br>
0_1800 : 269011<br>
2_1800 : 279397<br>
1_1800 : 440793<br>
2_60   : crashed<br>
0_900  : 504135<br>
2_900  : 533992<br>
0_60   : crashed<br>
3_3600 : 686641<br>
1_900  : 935043<br>
1_60   : crashed<br>
3_1800 : 1817668<br>
3_60   : stopped on January 10th 11.03am<br>
3_900  : stopped on January 10th 11.03am<br>
error margin 0.000001<br>
0_60   : stopped on January 10th 11.03am<br>
1_60   : crashed<br>
2_60   : crashed<br>

Three experiments were restarted with bigger error margins. They had crashed in earlier experiments. Due to the number of intervals, they encountered divide by zero errors. (Jaay for rounding and accuracy of computers)

For questions, feel free to contact Leoni Winschermann (L d o t Winschermann a t utwente d o t nl)

### Competitive ratio experiments
Experimental setup February 2024 <br>
Run with Python 3.8.8 <br>
We conducted experiments with the following parameters:

reps = 500<br>
instanceSizes = [400] <br>
timeSteps = [900]<br>
maxFlows = [shortest_augmenting_path]<br>
write = True<br>
randomSample = True <br>

The goal of the experiments was to put theoretical competitive ratios of two online algorithms into perspective. Here, an online problem is a problem where the existence, departure time, energy requirement and maximum charging power of a job becomes known upon its arrival. Therefore, we make scheduling decisions with the jobs currently available without knowledge of what other jobs may arrive after. <br>
The first considered online algorithm is Average Rate where upon arrival a job is scheduled at the average power needed to charge its energy requirement before departure. Assuming the 2-norm objective function, this algorithm has a theoretical competitive ratio of 8. <br>
The second considered online algorithm is Optimal Available where upon arrival the charging of all energy of already available jobs that has not yet been provided is reoptimized and the schedule updated accordingly. Assuming the 2-norm objective function, this algorithm has a theoretical competitive ratio of 4. <br>
The main file can be found under 'competitive_ratio_experiments/'. 

For questions, feel free to contact Leoni Winschermann (L d o t Winschermann a t utwente d o t nl)

### Carbon steering experiments
Experimental setup April 2024 <br>
Run with Python 3.11.7 <br>
Initial commit d4812cb <br>
We conducted experiments with the following parameters:

timeSteps = [900]<br>
maxFlows = [shortest_augmenting_path]

The used data stems from an electric bus charging hub. Therefore, an additional data processing file is included in carbon_steering_experiments/busses.py. The goal of the experiments was to validate a transform applied to incorporate carbon steering (similarly price steering) objectives in weighted combination with the original flattening objective into FOCS. Furthermore, the tradeoff between flatness and carbon intensity or the resulting schedule may be investigated based on the weights attached to either objective. <br>
The main file can be found under 'carbon_steering_experiments/'

For questions, feel free to contact Leoni Winschermann (L d o t Winschermann a t utwente d o t nl) or Leander van der Bijl (L d o t C d o t vanderBijl a t utwente d o t nl)

## References

[1] Leoni Winschermann, Marco E.T. Gerards, Antonios Antoniadis, Gerwin Hoogsteen, Johann Hurink. 2023. Relating Electric Vehicle Charging to Speed Scaling with Job-Specific Speed Limits. https://arxiv.org/abs/2309.06174. [currently under peer review] <br>
[2] Aric A. Hagberg, Daniel A. Schult, and Pieter J. Swart. 2008. Exploring Network Structure, Dynamics, and Function using NetworkX. In Proceedings of the 7th Python in Science Conference, Gaël Varoquaux, Travis Vaught, and Jarrod Millman (Eds.). Pasadena, CA USA, 11 – 15. <br>
[3] Gurobi Optimization, LLC. 2023. Gurobi Optimizer Reference Manual. https://www.gurobi.com <br>

# Related papers

Leoni Winschermann, Marco E.T. Gerards, Antonios Antoniadis, Gerwin Hoogsteen, Johann Hurink. 2023. Relating Electric Vehicle Charging to Speed Scaling with Job-Specific Speed Limits. https://arxiv.org/abs/2309.06174. [currently under peer review] <br>
Leoni Winschermann, Marco E.T. Gerards, Johann Hurink. 2024. Improving the Optimization in Model Predictive Controllers: Scheduling Large Groups of Electric Vehicles. https://arxiv.org/abs/2403.16622. [currently under peer review] <br>
