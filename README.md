[![Python 3.7](https://img.shields.io/badge/python-3.7-green.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

[![Contact Info](https://img.shields.io/badge/Contact%20Info-fangfeili.fanny@gmail.com-orange.svg)]()
[![Update Info](https://img.shields.io/badge/Update%20Info-Feb20,%202020-orange.svg)]()


## PyFitMut

### 1. What is PyFitMut?

PyFitMut is a Python-based tool that can identify spontaneous adaptive mutations for initially isogenic evolving population, as well as estimate the fitness effect and establishment time of those adaptive mutations. The detailed theory and algorithm of PyFitMut is introduced in reference: [S. F. Levy, et al. Quantitative Evolutionary Dynamics Using High-resolution Lineage Tracking. Nature, 519: 181-186 (2015)](https://www.nature.com/articles/nature14279). If you use this software, please reference: [???](???). PyFitMut is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

It currently has two main functions:
* `evomut_simulator.py` performs simulations of competitve pooled growth of a isogenic population with spotaneous mutations.
* `pyfitmut.py` calculates the fitness and establishment time of all adaptive mutations from read-count time-series data.
    
A walk-through is included as the jupyter notebook [here](https://github.com/FangfeiLi05/PyFitSeq/blob/master/PyFitMut_Walk_Through.ipynb).



### 2. How to install PyFitMut?
* Python 3 is required. This version has been tested on a MacBook Pro (3.1 GHz Intel Core i5), with Python 3.7.4.
* Clone this repository by running `git clone https://github.com/FangfeiLi05/PyFitMut.git` in terminal.
* `cd` to the root directory of the project (the folder containing `README.md`).
* Install dependencies by running `pip install -r requirements.txt` in terminal.



### 3. How to use PyFitMut?

#### 3.1. Evolution Simulation
`evomut_simulator.py` models competative pooled growth of a isogenic population with spotaneous mutations. This simulation can be made to include an arbitry probability distribution for finess effect of mutations，and include sources of noise, including growth noise, noise from cell transfers, DNA extraction, PCR, and sequencing.

##### OPTIONS
* `--input` or `-i`: a .csv file, with
  + 1st column of .csv: initial cell number of each genotype at generation 0, [n1, n2, ...]
  + 2nd column of .csv: bin edges of the arbitrary probability distribution for fitness effect of mutations
  + 3rd column of .csv: counts frequency in each bin of the arbitrary probability distribution for fitness effect of mutations
* `--t_seq` or `-t`: time-points evaluated in number of generations (`format: 0 t1 t2 ...`)
* `--read_num_average_seq` or `-r`: average number of reads per genotype for each time-point (`format: 0 r1 r2 ...`)
* `--noise_option` or `-n`: which types of noise to include in the simulation, default is all sources of noise (`default: growth bottleneck_transfer DNA_extraction PCR sequencing`)
* `--dna_copies` or `-d`: average genome copy number per genotype used as template in PCR (`default: 500`)
* `--pcr_cycles` or `-p`: number of cycles of PCR (`default: 25`)
* `--maximum_mutation_number` or `-m`: number of maximum mutations allowed for each lineage (`default: 1`)
* `--mutation_rate` or `-u`: total beneficial mutation rate (`default: 1e-5`)
* `--output_filename` or `-o`: prefix of output .csv files (`default: output`)

##### OUTPUTS
* `output_filename_EvoSimulation_Read_Number.csv`: read number per genotype for each time-point
* `output_filename_EvoSimulation_Other_Info.csv`: information for mutations (include fitness effect, occured time, establishment time, lineage label, mutation label within a lineage), mean fitness for each time-point, a record of all inputs(including t_seq, read_num_average_seq, noise_option, dna_copies, maximum_mutation_number, mutation_rate)




##### For Help
```
python mutevo_simulator.py --help
```

##### Examples
```
python mutevo_simulator.py -i simu_input_exp.csv -t 0 9 18 27 36 45 54 63 72 81 90 99 108 -r 50 50 50 50 50 50 50 50 50 50 50 50 50
python mutevo_simulator.py -i simu_input_exp.csv -t 0 9 18 27 36 45 54 63 72 81 90 99 108 -r 75 75 75 75 75 75 75 75 75 75 75 75 50 -n DNA_extraction PCR sequencing -d 300 -p 27 -u 2e-5 -m 2 -o output
```      



#### 3.2. Fitness Estimation
`pyfitmut.py` identifies adaptive mutations and estimates their fitness effect as well as establishment time from read-count time-series data.


##### OPTIONS
* `--input` or `-i`: a .csv file, with each column being the read number per genotype at each sequenced time-point
* `--t_seq` or `-t`: sequenced time-points in number of generations (`format: 0 t1 t2 ...`)
* `--output_filename` or `-o`: prefix of output .csv files (`default: output`)

##### OUTPUTS
* `output_filename_FitSeq_Result.csv`: a .csv file, with
  + 1st column of .csv: estimated fitness of each lineage (0 if there is no mutation), [x1, x2, ...]
  + 2nd column of .csv: estimated estabblishment time of each lineage (0 if there is no mutation), [tau1, tau2, ...]
  + 3rd column of .csv: log likelihood value of each lineage (0 if there is no mutation), [f1, f2, ...]
  + 4rd column of .csv: estimated mean fitness per sequenced time-point, [x_mean(0), x_mean(t1), ...]

##### For Help
```
python pyfitmut.py --help
```  

##### Examples
```
python pyfitmut.py -i output_EvoSimulation_Read_Number.csv -t 0 9 18 27 36 45 54 63 72 81 90 99 108 -n 9e8 9e8 9e8 9e8 9e8 9e8 9e8 9e8 9e8 9e8 9e8 9e8 9e8 -u 2e-5
python pyfitmut.py -i output_EvoSimulation_Read_Number2.csv -t 0 2 6 8 -m 12 -k 2 -g 3 -f w -o output
``` 




