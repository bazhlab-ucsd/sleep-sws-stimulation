
# Stimulation Augments Spike Sequence Replay and Memory Consolidation during Slow-Wave Sleep

This project reproduces the results of [Stimulation Augments Spike Sequence Replay and Memory Consolidation during Slow-Wave Sleep](https://www.jneurosci.org/content/40/4/811).

## Usage

The code is implemented in C++ and uses OpenMP multi-core parallelism to optimize the code's performance.
To generate the network connectivity file, navigate to the root directory of the project and run the following command: 

`make network`

To run the simulation, use the following command:

`make run` 

The code generates several output files in the out folder that illustrate the electrical activity in the brain under different levels of neuromodulators. The output files include the membrane voltage of cortical neurons (time_cx) and other neuron types. You can vary the acetylcholine, histamine, GABA, and monoamine levels by modifying the parameters in the params.txt file. The network connectivity is given in the `network.cfg` file.

The `network.cfg` file includes the connectivity information of different brain regions and neuron types. The `params.txt` file contains several parameters, such as the levels of neuromodulators, the time step, and the simulation duration. You can modify these parameters to customize the simulation.

The output files are saved in the out folder with a timestamp in the filename. Each file contains the activity of different neurons at different times during the simulation. The membrane voltage of cortical neurons is saved in the time_cx file. Other neuron types are saved in separate files with their respective names.


### About the model

The computer model used in the study considers computational models implementing effects of neuromodulators to simulate transitions between awake and SWS sleep, and synaptic plasticity to allow the change of synaptic connections due to the training in awake or replay during sleep. 

The default parameters are set to replicate the results for figure 2 in the following publication.


### Please refer to following publication for more details:

Stimulation Augments Spike Sequence Replay and Memory Consolidation during Slow-Wave Sleep.<br>
Yina Wei, Giri P. Krishnan, Lisa Marshall, Thomas Martinetz and Maxim Bazhenov. Journal of Neuroscience 22 January 2020, 40 (4) 811-824; [DOI: https://doi.org/10.1523/JNEUROSCI.1427-19.2019 ](https://doi.org/10.1523/JNEUROSCI.1427-19.2019)
