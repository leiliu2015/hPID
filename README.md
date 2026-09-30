## hPID

Two disparate mechanisms, namely *specific* button model versus *nonspecific* button model, have been proposed to explain robust somatic homolog pairing in *Drosophila*; yet, how homolog pairing shapes 3D genome remains poorly understood. Here, we performed Brownian dynamics simulations of these two models by considering the following papers:
- M. B. Child VI, J. R. Bateman, A. Jahangiri, A. Reimer, N. C. Lammers, N. Sabouni, D. Villamarin, G. C. McKenzie-Smith, J. E. Johnson, D. Jost, and H. G. Garcia, [Live imaging and biophysical modeling support a button-based mechanism of somatic homolog pairing in *Drosophila*](https://elifesciences.org/articles/64412), eLife 10, e64412 (2021).
- W. F. Marshall, and J. C. Fung, [Modeling homologous chromosome recognition via nonspecific interactions](https://www.pnas.org/doi/10.1073/pnas.2317373121), PNAS 121, e2317373121 (2024).

Based on the results, a *bubble-chain* model is proposed which clarifies the geometrical origin of homolog pairing-induced contact domains (***hPID****s*) in the *Drosophila* diploid genome.

### System Requirements

Codes written in Python and [Gnuplot](https://gnuplot.sourceforge.net/) to reproduce most results in our recent [work]() are archived in this repository, which have been tested on ubuntu 18.04/20.04 LTS. We recommend [Anaconda](https://www.anaconda.com/) to manage the Python environment (*>=3.7*) and packages, including [Numpy](https://numpy.org) and [Numba](numba.pydata.org). The latter is essential to shorten the running time of simulations. 

### File Description
- codes/
  - [example30.seq](codes/example30.seq) (A TXT file containing two lines, each specifying a button sequence for a pair of homologous chromatin fibers, where numbers 0 and 1 represent nonadhesive and button loci, respectively)
  - [example10.seq](codes/example30.seq) (A similar TXT file containing two button sequences with a reduced linear density, which relates to Fig. 6 and Fig. S6)
  - [pid.nb.py](codes/pid.nb.py) (A Python script to simulate *nonspecific* button model)
  - [pid.sb.py](codes/pid.sb.py) (A Python script to simulate *specific* button model)
  - [pid.bb.py](codes/pid.bb.py) (A Python script to simulate *bubble-chain* model)
  - [pid.xb.xyz2h5.py](codes/pid.xb.xyz2h5.py) (A Python script to change simulation trajectories into HDF5 format)
  - [pid.xb.h5px.dist.py](codes/pid.xb.h5px.dist.py) (A Python script to calculate the contact map based on trajectories)
  - [mdRandomSeeds](codes/mdRandomSeeds) (A TXT file of random number seeds)
  - [reds.pal](codes/reds.pal) (A Gnuplot color palette file)
- nonspecific-button/
  - [tutorial-nb.sh](nonspecific-button/tutorial-nb.sh) (A BASH script to perform all the modeling and analyses in this directory)
  - [tutorial-nb.gnu](nonspecific-button/tutorial-nb.gnu) (A Gnuplot script to visualize the modeling results)
  - backup/ (A backup of the simulation results for reference)
- other directories/
  - Similar to [nonspecific-button/](nonspecific-button/), each directory contains one BASH script and one Gnuplot script for one different model. 
- [clearAll.sh](clearAll.sh) (A BASH script to delete all outputs) 

### User Guide

Taking nonspecific button model as an example, first, open your terminal and run the BASH script 
```
$ cd ./nonspecific-button
$ bash ./tutorial-nb.sh
```
It takes less than 20 minutes to finish on our desktop with a Intel® Core™ i5-9500 CPU. Next, one can plot the results by using the Gnuplot script
```
$ gnuplot -persist tutorial-nb.gnu
```
Each simulation is composed of two stages, and the button-based mechanisms are switched on *only* in the second stage. The Gnuplot script shows the time dependence of the end-to-end distance of the polymer chains (*R<sub>ee</sub>*), the fraction of correctly paired buttons (*f<sub>pair</sub>*), the final configurations of the system, and the time-averaged contact probability map (*p<sub>ij</sub>*). To plot results in the backup folder, set the variable `fd` in the first line of [tutorial-nb.gnu](nonspecific-button/tutorial-nb.gnu) to `'backup/'`. 

The other two models can be simulated in a similar manner. The approximate running time for each tutorial and their related figure indices in our [preprint]() are listed in the following table. *Note* that the total simulation time of these tutorials, which can be specified by changing the value of the variable `numberOfSteps` in the modeling scripts, has been reduced to 1/50 of that in the paper. 

| Directory | Running Time (min) | Figure |
| --------- | ------------------ | -------|
|nonspecific-button | 20 | Fig. 3|
|specific-button | 20 | Fig. 4 and S3|
|bubble-chain | 15 | Fig. S4|

Last, all the output files can be deleted by typing `$ bash ./clearAll.sh` at the root of this repository. For any inquiry about these codes, feel free to contact Lei Liu (leiliu2015@163.com)

