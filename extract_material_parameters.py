import numpy as np
from scipy.optimize import curve_fit

# Class to store data from a single experiment
class experiment():
    def __init__(self,data,pressure=5):
        # Initialize the experiment with the data and pressure
        self.t = data[:,0]
        self.L = data[:,1]

        # Remove NaN values
        nan_indicies = np.isnan(self.L)
        self.L = self.L[~nan_indicies]
        self.t = self.t[~nan_indicies]

        # Set the initial length and remove the initial length from the data
        self.L0 = self.L[0]
        self.L = self.L - self.L0

        # Set the maximum length
        self.Lmax = np.max(self.L)

        # Set the pressure
        self.pressure = pressure

# Import data from a csv file and return a list of experiments
def import_data(file_name,pressure=5):
    data = np.genfromtxt(file_name,skip_header=1,delimiter=',',filling_values=np.nan)
    experiment_list = []
    for i in range(0,data.shape[1],2):
        experiment_list.append(experiment(data[:,i:i+2],pressure=pressure))
    return experiment_list

# Import data from csv file with a list of pressures in the header
def import_data_laplace(file_name):
    data = np.genfromtxt(file_name,skip_header=1,delimiter=',',filling_values=np.nan)
    with open(file_name) as f:
        first_line = f.readline().strip('\n')
    first_line = first_line.split(',')[1:]
    pressures = []
    for entry in first_line:
        p,_ = entry.split(' ')
        pressures.append(float(p))

    n = data.shape[1]-1
    experiment_list = []
    for i in range(1,n+1):
        experiment_list.append(experiment(data[:,[0,i]],pressure=pressures[i-1]))
    return experiment_list

# Below is code for fitting expeirmental Lp(t) data to a given model and extracting material parameters
# Define models for solid and fluid like behavior. Models constrained as we assume initial length is zero
def constrained_solid_model(t,A,D):
    return A*(1-np.exp(-D*t))

def constrained_fluid_model(t,B):
    return B*t

def constrained_general_model(t,A,B,D):
    return constrained_fluid_model(t,B) + constrained_solid_model(t,A,D)

def weighted_mean_and_standard_deviation(x,x_std):
    x_mean = np.sum(x/x_std**2)/np.sum(1/x_std**2)
    x_mean_std = np.sqrt(1/np.sum(1/x_std**2))
    return x_mean,x_mean_std

def error_propogation_through_mean(x,x_std):
    x_mean = np.mean(x)
    x_mean_std = np.sqrt(np.sum(x_std**2))/len(x_std)
    return x_mean,x_mean_std

def mean_and_standard_error(x,x_std):
    x_mean = np.mean(x)
    x_sample_std = np.std(x)
    x_sem = x_sample_std/np.sqrt(len(x))
    return x_mean,x_sem

def sample_mean_and_standard_deviation(x,x_std):
    x = x[x!=0]
    x_std = x_std[x_std!=0]
    return mean_and_standard_error(x,x_std)

# Fit an experiment to a given model:
def fit_experiment(experiment,model,guess,Rp=1e-6):
    # This accounts for the 0.2 um uncertainity in measurements of Lp and results in higher uncertainity in final results
    measurement_error = np.ones(experiment.t.shape)*(1/Rp*.2e-6)
    pout,pcov = curve_fit(model,experiment.t,experiment.L,p0=guess,maxfev=8000,sigma=measurement_error,absolute_sigma=True)
    pvar = np.diagonal(pcov)
    return pout, pvar

# Fit a list of experiments to a given model:
def fit_experiment_list(experiment_list,model,guess):
    popt_list = []
    pvar_list = []
    for experiment in experiment_list:
        popt,pvar = fit_experiment(experiment,model,guess)
        popt_list.append(popt)
        pvar_list.append(pvar)
    popt_list = np.squeeze(np.array(popt_list))
    pvar_list = np.squeeze(np.array(pvar_list))
    return popt_list,pvar_list

# This function takes a set of fitting parameters and their variances and groups them by sample
# The parameters are averaged over the samples and the variance is propagated to find the standard deviation
def group_fitting_parameters_by_sample(popt_list,pvar_list,groups):
    # We first build a mask for each group to enable propertly averaging over the sample nucleoli
    num_groups = np.max(groups)+1
    n_samples = len(groups)
    mask = np.array([(n_samples*[i]) == groups for i in range(num_groups)]).astype(int)
     # Apply the mask to the parameters
    # Average over each nucleoli sample and propagate the variance to find the standard deviation or standard error
    # When taking the mean of samples with individual variances, the standard deviation of the mean is given by the square root of the sum of the variances divided by the number of samples (propogation of uncertainty) rather than the traditional standard error of the mean formula
    if popt_list.ndim == 1:
        popt_list = mask * popt_list
        pvar_list = mask * pvar_list
        popt_list = popt_list[...,np.newaxis]
        pvar_list = pvar_list[...,np.newaxis]
    else:
        popt_list = mask[...,np.newaxis] * popt_list
        pvar_list = mask[...,np.newaxis] * pvar_list
    mean_popt_list = []
    mean_popt_std_list = []
    for i in range(num_groups):
        val, var = [], []
        for j in range(popt_list.shape[2]):
            mean_val, mean_var = sample_mean_and_standard_deviation(popt_list[i,:,j],pvar_list[i,:,j])
            val.append(mean_val)
            var.append(mean_var)
        mean_popt_list.append(val)
        mean_popt_std_list.append(var)
    mean_popt = np.squeeze(np.array(mean_popt_list))
    mean_popt_std = np.squeeze(np.array(mean_popt_std_list))
    return mean_popt,mean_popt_std

# Compute the viscosity and its standard deviation from fitting parameters generated by a constrained fluid model
# The input parameters are the fitting parameters, their standard deviations, the radius of the pipette, the radius of the sample, the applied pressure, the surface tension, and the standard deviation of the surface tension
# Optionally, we can supply a mask to exclude certain samples from the calculation
def calculate_fluid_properties(popt_list,pstd_list,rp,rc,pressure,gamma,gamma_std,mask=None):
    # Calculate the Laplace pressure
    laplace_pressure = 2*gamma*(1/rp-1/rc)
    delta_pressure = pressure - laplace_pressure
    # Calculate STD of delta pressure assuming 1Pa error in pressure and following uncertainty propogation rules
    pressure_std = np.sqrt((1.)**2*np.ones_like(delta_pressure)+(gamma_std*2*(1/rp-1/rc))**2)
    # Calculate the viscosity
    viscosity = delta_pressure/popt_list
    # Calculate the standard deviation of the viscosity
    viscosity_std = np.sqrt(viscosity**2*((pressure_std/delta_pressure)**2+(pstd_list/popt_list)**2))
    # Apply the mask to the viscosity and viscosity standard deviation
    if mask is not None:
        viscosity = viscosity[mask]
        viscosity_std = viscosity_std[mask]
    # Calculate the mean viscosity and its standard deviation
    mean_viscosity,mean_viscosity_std = sample_mean_and_standard_deviation(viscosity,viscosity_std)
    return mean_viscosity,mean_viscosity_std

def calculate_solid_properties(popt_list,pstd_list,rp,pressure,gamma,gamma_std,mask=None):
    # Calculate the Laplace pressure, noting that no DFC samples have a measured rc value
    laplace_pressure = 2*gamma*(1/rp)
    delta_pressure = pressure - laplace_pressure
    # Calculate STD of delta pressure assuming 1Pa error in pressure and following uncertainty propogation rules
    pressure_std = np.sqrt((1.)**2*np.ones_like(delta_pressure)+(gamma_std*2*(1/rp))**2)
    # Calculate the elastic modulus
    elastic_modulus = delta_pressure/popt_list[:,0]
    # Calculate the standard deviation of the elastic modulus
    elastic_modulus_std = np.sqrt(elastic_modulus**2*((pstd_list[:,0]/popt_list[:,0])**2+(pressure_std/delta_pressure)**2))
    # Calculate the viscosity
    viscosity = elastic_modulus/popt_list[:,1]
    # Calculate the standard deviation of the viscosity following propogation of uncertainity
    viscosity_std = np.sqrt(viscosity**2*((elastic_modulus_std/elastic_modulus)**2+(pstd_list[:,1]/popt_list[:,1])**2))
    # Apply the mask to the properties and their standard deviations
    if mask is not None:
        viscosity = viscosity[mask]
        viscosity_std = viscosity_std[mask]
        elastic_modulus = elastic_modulus[mask]
        elastic_modulus_std = elastic_modulus_std[mask]
    # Calculate the mean properties and their standard deviation
    mean_viscosity,mean_viscosity_std = sample_mean_and_standard_deviation(viscosity,viscosity_std)
    mean_elastic_modulus,mean_elastic_modulus_std = sample_mean_and_standard_deviation(elastic_modulus,elastic_modulus_std)
    return mean_elastic_modulus,mean_elastic_modulus_std,mean_viscosity,mean_viscosity_std

# Below is code for fitting a laplace pressure experiment and extracting the surface tension
def linear_response(x,A,B):
    return A*x + B

def fit_laplace_pressure_experiment(file,full_return=False,Rp=1):
    exp_list = import_data_laplace(file)
    pressures = []
    popt_list = []
    pvar_list = []
    # Fit each experiment to a constrained fluid model and store the parameters
    for exp in exp_list:
        popt,pvar = fit_experiment(exp,constrained_fluid_model,[1],Rp=Rp)
        pressures.append(exp.pressure)
        popt_list.append(popt)
        pvar_list.append(pvar)
    popt_list = np.squeeze(np.array(popt_list))
    pvar_list = np.squeeze(np.array(pvar_list))
    p_std_list = np.sqrt(pvar_list)

    # Perform linear regression on the parameters, from which the surface tension can be extracted as the y-intercept
    popt,pcov = curve_fit(linear_response,pressures,popt_list,sigma=p_std_list,absolute_sigma=False)
    pstd = np.sqrt(np.diagonal(pcov))
    # Compute the y-intercept and its standard deviation
    p_intercept = - popt[1]/popt[0]
    p_std = np.sqrt(p_intercept**2*((pstd[0]/popt[0])**2 + (pstd[1]/popt[1])**2 + 2*pcov[0,1]/(popt[0]*popt[1])))
    return p_intercept,p_std

def fit_surface_tension(file_list,Rp,Rc):
    p_intercept_list = []
    p_std_list = []
    # Find the y-intercept of the linear regression on the Laplace pressures for each experiment
    for file,Rp_val in zip(file_list,Rp):
        p_intercept,p_std = fit_laplace_pressure_experiment(file,Rp=Rp_val)
        p_intercept_list.append(p_intercept)
        p_std_list.append(p_std)
    p_intercept = np.array(p_intercept_list)
    p_var = np.array(p_std_list)**2

    # Compute the surface tension and its standard deviation
    pre_factor = 0.5*np.divide(np.ones_like(Rc),np.divide(np.ones_like(Rp),Rp) - np.divide(np.ones_like(Rc),Rc))
    gamma = pre_factor*p_intercept
    gamma_variance = pre_factor**2 * (p_var)
    gamma_mean,gamma_std = sample_mean_and_standard_deviation(gamma,np.sqrt(gamma_variance))
    return gamma_mean,gamma_std

if __name__ == "__main__":
    # Define GC Laplace Experiment Data:
    gc_laplace_experiment_file_list = [f'laplace_pressure_gc/gc_{i}.csv' for i in range(1,5)]
    gc_laplace_Rp = np.array([1,1,1,.75])*1e-6
    gc_laplace_Rc = np.array([5.3,3.5,5.1,3.6])*1e-6
    gc_gamma, gc_gamma_std = fit_surface_tension(gc_laplace_experiment_file_list,gc_laplace_Rp,gc_laplace_Rc)

    # Define GC data and magic numbers
    gc_experiment_list = import_data("compiled_data/compiled_gc.csv")
    gc_groups = np.array([0,1,1,1,2,2,2,3,4,4,4,5,6,7])
    gc_rp = np.array([1.00E-06,1.00E-06,1.00E-06,1.00E-06,6.50E-07,7.50E-07,7.00E-07,8.00E-07])
    gc_rc = np.array([0.0000053,3.50E-06,5.10E-06,0.00000475,5.00E-06,0.0000036,0.0000042,0.0000052])
    # This mask excludes all samples that are over 60 minutes old
    gc_mask = np.array([True,True,False,False,True,True,False,True])
    gc_pressure = 5

    # Fit the GC data to a constrained fluid model and extract the fluid properties
    gc_popt_list, gc_pvar_list = fit_experiment_list(gc_experiment_list,constrained_fluid_model,[.1])
    gc_popt_list, gc_pstd_list = group_fitting_parameters_by_sample(gc_popt_list,gc_pvar_list,gc_groups)
    gc_viscosity,gc_viscosity_std = calculate_fluid_properties(gc_popt_list,gc_pstd_list,gc_rp,gc_rc,gc_pressure,gc_gamma,gc_gamma_std)
    young_gc_viscosity, young_gc_viscosity_std = calculate_fluid_properties(gc_popt_list,gc_pstd_list,gc_rp,gc_rc,gc_pressure,gc_gamma,gc_gamma_std,mask=gc_mask)

    # Define DFC data and magic numbers
    dfc_experiment_list = import_data("compiled_data/compiled_dfc.csv")
    dfc_groups = np.array([0,1,1,2,2,3,3,4,5,6,7,8,9])
    dfc_rp = np.array([9.00E-07,4.50E-07,4.50E-07,5.00E-07,8.50E-07,8.50E-07,7.00E-07,6.00E-07,5.00E-07,8.00E-07])
    dfc_pressure = 20
    dfc_gamma = gc_gamma
    dfc_gamma_std = gc_gamma_std

    # Fit the DFC data to a constrained solid model and extract the solid properties
    dfc_popt_list, dfc_pvar_list = fit_experiment_list(dfc_experiment_list,constrained_solid_model,[1,.01])
    dfc_popt_list, dfc_pstd_list = group_fitting_parameters_by_sample(dfc_popt_list,dfc_pvar_list,dfc_groups)
    dfc_elastic_modulus, dfc_elastic_modulus_std, dfc_viscosity, dfc_viscosity_std = calculate_solid_properties(dfc_popt_list,dfc_pstd_list,dfc_rp,dfc_pressure,dfc_gamma,dfc_gamma_std)

    # Define RNase treated DFC Laplace Experiment Data:
    rnase_files = np.asarray([f'laplace_pressure_rnase/rnase_{i}.csv' for i in range(1,5)])
    rnase_Rp = np.array([0.75,.75,.6,.6])*1e-6
    rnase_Rc = np.array([1.75,1.8,1.5,2.3])*1e-6
    # rnase_files = rnase_files[[1,2,3]]
    # rnase_Rp = rnase_Rp[[1,2,3]]
    # rnase_Rc = rnase_Rc[[1,2,3]]
    rnase_gamma, rnase_gamma_std = fit_surface_tension(rnase_files,rnase_Rp,rnase_Rc)

    # Define RNase treated DFC data and magic numbers
    rnase_experiment_list = import_data("compiled_data/compiled_rnase.csv")
    rnase_groups = np.array([0,1,2,2,3,4,4,4,4])
    rnase_rp = np.array([7.50E-07,7.50E-07,7.50E-07,6.00E-07,6.00E-07])
    rnase_rc = np.array([1.80E-06,1.90E-06,1.75E-06,2.30E-06,1.50E-06])
    rnase_pressure = 20

    # Fit the RNase treated DFC data to a constrained fluid model and extract the fluid properties
    rnase_popt_list, rnase_pvar_list = fit_experiment_list(rnase_experiment_list,constrained_fluid_model,[.1])
    rnase_popt_list, rnase_pstd_list = group_fitting_parameters_by_sample(rnase_popt_list,rnase_pvar_list,rnase_groups)
    rnase_viscosity,rnase_viscosity_std = calculate_fluid_properties(rnase_popt_list,rnase_pstd_list,rnase_rp,rnase_rc,rnase_pressure,rnase_gamma,rnase_gamma_std)

    # Print the results
    print("GC Material Parameters:")
    print(f"GC Viscosity: {gc_viscosity:.2e} ± {gc_viscosity_std:.2e} Pa·s")
    print(f"Young GC Viscosity: {young_gc_viscosity:.2e} ± {young_gc_viscosity_std:.2e} Pa·s")
    print(f"GC Surface Tension: {gc_gamma:.2e} ± {gc_gamma_std:.2e} N/m")
    print("DFC Material Parameters:")
    print(f"DFC Elastic Modulus: {dfc_elastic_modulus:.2e} ± {dfc_elastic_modulus_std:.2e} Pa")
    print(f"DFC Viscosity: {dfc_viscosity:.2e} ± {dfc_viscosity_std:.2e} Pa·s")
    print("RNase Material Parameters:")
    print(f"RNase Viscosity: {rnase_viscosity:.2e} ± {rnase_viscosity_std:.2e} Pa·s")
    print(f"RNase Surface Tension: {rnase_gamma:.2e} ± {rnase_gamma_std:.2e} N/m")