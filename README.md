# cheng_et_al_2025_data
This repository contains the data presented in ``Micropipette aspiration reveals differential RNA-dependent viscoelasticity of nucleolar subcompartments'' by Holly H. Cheng, James V. Roggeveen, Huan Wang, Howard A. Stone, Zheng Shi, Clifford P. Brangwynne. Accompanying this data is a python script `extract_material_parameters.py` that can be used to extract the material parameters from the data and to compute the errors assocaited with these parameters following standard error propagation rules.

## Data
The data is present in three folders:
- `compiled_data`
- `laplace_pressure_gc`
- `laplace_pressure_rnase`

In `compiled_data` there are three csv files, which contain the data for the GC, DFC, and RNase treated DFC experiments. Each of these files consists of 2N columns, where N is the number of discrete experimental segments. The columns repeat as time then Lp(t). The GC segments were all taken at a pressure of 5 Pa while the DFC and RNase treated DFC segments were taken at 20 Pa.

The `laplace_pressure_gc` and `laplace_pressure_rnase` folders contain the data for the Laplace pressure experiments for the GC and RNase treated DFC experiments, respectively. Each of these folders contains a set of csv files, each of which contain the data for a single experiment. These files have the form time then Lp(t) for several pressures. The pressure for each column is indicated in the header.

## Analysis
The `extract_material_parameters.py` script can be used to extract the material parameters from the data and to compute the errors assocaited with these parameters following standard error propagation rules. The script can be run by executing the command `python extract_material_parameters.py` from the root directory. The parameters will print to the terminal. A copy of the most recent output can be found in `material_parameters.txt`.

## Requirements
The script was written in python 3.13 and has the following dependencies:
- numpy
- scipy