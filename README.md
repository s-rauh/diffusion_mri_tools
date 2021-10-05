# Diffusion MRI Tools

Diffusion MRI Tools contains Matlab code to analyze and process diffusion MRI data. 

Three types of fitting are currently implemented: 
* IVIM fit
* DTI fit
* IVIM-DTI fit

### IVIM fit
The IVIM model estimates the tissue diffusion coefficient D, the perfusion fraction f and the pseudo-diffusion coefficient D*. Three fit methods are supported: free, two-step approach or segmented/simplified IVIM. 

### DTI fit
The diffusion tensor is estimated from the provided data. A contrained and an unconstrained fit are implemented. The constrained fit uses a Cholesky decomposition to ensure the diffusion tensor is positive semi-definite.  

### IVIM-DTI fit
A combined IVIM-DTI fit is performed. The following methods are currently implemented: Free fit, two-step fit, segmented/simplified IVIM-DTI fit, IVIM-corrected DTI fit. 

## Analysis tools
The function calc_dti_parameters calculated the MD, FA, eigenvalues and eigenvectors from a diffusion tensor. 
