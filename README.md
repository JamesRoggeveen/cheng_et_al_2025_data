# cheng_et_al_2025_data
This repository contains the data presented in ``Micropipette aspiration reveals differential RNA-dependent viscoelasticity of nucleolar subcompartments'' by Holly H. Cheng, James V. Roggeveen, Huan Wang, Howard A. Stone, Zheng Shi, Clifford P. Brangwynne. Accompanying this data is a python script `extract_material_parameters.py` that can be used to extract the material parameters from the data and to compute the errors associated with these parameters following standard error propagation rules.

## Data
The data is present in three folders:
- `compiled_data`
- `laplace_pressure_gc`
- `laplace_pressure_rnase`

In `compiled_data` there are three csv files, which contain the data for the GC, DFC, and RNase treated DFC experiments. Each of these files consists of 2N columns, where N is the number of discrete experimental segments. The columns repeat as time then Lp(t). The GC segments were all taken at a pressure of 5 Pa while the DFC and RNase treated DFC segments were taken at 20 Pa.

The `laplace_pressure_gc` and `laplace_pressure_rnase` folders contain the data for the Laplace pressure experiments for the GC and RNase treated DFC experiments, respectively. Each of these folders contains a set of csv files, each of which contain the data for a single experiment. These files have the form time then Lp(t) for several pressures. The pressure for each column is indicated in the header.

## Analysis
The `extract_material_parameters.py` script can be used to extract the material parameters from the data and to compute the errors associated with these parameters following standard error propagation rules. The script can be run by executing the command `python extract_material_parameters.py` from the root directory. The parameters will print to the terminal.

Within the script, the function `sample_mean_and_standard_deviation` can be modified depending on the desired method of averaging the data and computing the resulting uncertainties by changing the function call prior to the return statement. Three options are available:
- `mean_and_standard_error`: This option computes the mean and standard error of the data. This option does not account for the uncertainty in the measurement process.
- `error_propogation_through_mean`: This option computes the mean and the standard deviation of the data assuming that the data is normally distributed and that the uncertainty is due to the measurement process. This option does not account for the variation between samples.
- `weighted_mean_and_standard_deviation`: This option computes the mean weighted by the inverse of the variance of the data and the corresponding standard deviation. This is the minimum variance estimator of the mean.

The file `material_parameters_error_propagation.txt` contains the most recent material parameters and their errors computed using the `error_propogation_through_mean` function. The file `material_parameters_standard_error.txt` contains the most recent material parameters and their standard errors computed using the `mean_and_standard_error` function. In our paper, we report the error as the larger of the two error estimates.

## Requirements
The script was written in python 3.13 and has the following dependencies:
- numpy
- scipy