%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to perform an DTI fit to data. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input values: Diffusion data, b-values, diffusion directions
% Options: 
%   - mask:             Default 1
%                       Background is masked using a threshold cut-off 
%                       value 
%
%   - constrained:      Default 1
%                       Constrained: fit the lower triangular matrix of a
%                       Cholesky decomposition. The fit boundaries differ
%                       from the unconstrained fit. After the tensor fit,
%                       the diffusion tensor needs to be calculated from the
%                       lower triangular matrix. 
%                       Unconstrained: fit the 6 diffusion tensor elements
%                       immidiately. 
%
%   - normalized:       Default 1
%                       Normalize data to min(bval)-signal
%
% Output:
% Structure dti_fit containing the following fields:
%   - S0 signal (fitted). In case of normalization S0 will be close to 1
%     for all voxels. 
%   - diffusion tensor D
%     vector containing the 6 elements of the diffusion tensor. 
%     Order of the tensor elements is:
%     Dxx, Dyy, Dzz, Dxy, Dxz, Dyz. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dti_fit = fit_dti(data, bval, diffdir, options)

arguments
    data
    bval
    diffdir
    options.mask {mustBeNumericOrLogical} = 1
    options.constrained {mustBeNumericOrLogical} = 1
    options.normalize {mustBeNumericOrLogical} = 1
end

%rearrange and reshape data to a bval x n array
[data, sz] = reshape_diffdata_for_fit(data, bval);

if size(bval,1)==1
    bval = bval.';
end
%% optional: mask background voxels
%define datafit
if options.mask
    [datafit, sels] = mask_diffdata(data, bval);
else
    sels = true(1,size(data,2));
    datafit = data;
end

%% optional: normalize data
if options.normalize
    datafit = norm_diffdata(datafit, bval);
    %set initial guess
    x0 = [1 1 1 1 0 0 0];
else
    x0 = [mean(datafit(bval==min(bval),:), 'all') 1 1 1 0 0 0];
end

%% calculate b-matrix
bmat = calc_bmat(bval, diffdir);

%% set fit options
if options.constrained
    %Fit Cholesky lower triangular matrix of diffusion tensor
    lb = [0   -9 -9 -9 -9 -9 -9];
    ub = [inf  9  9  9  9  9  9];
else
    lb = [0    0  0  0  0  0  0];
    ub = [inf  3  3  3  3  3  3];
end

if options.normalize
    ub(1) = 10;
end
fitoptions = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt', ...
    'FunctionTolerance', 1e-10, 'MaxIterations', 1000, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
    'Display', 'off');

%initialize diffusion tensor and S0
tensor = zeros(6, size(datafit,2));
S0_fit = zeros(1, size(datafit,2));

%% perform fit
fprintf('Perform voxel-wise DTI fit...\n');
tic; 
for v = 1:size(datafit,2)
    if options.constrained
        t = lsqcurvefit(@dtifun_constr, x0, bmat, datafit(:,v), lb, ub, fitoptions);
        % Calculate tensor elements from Cholesky lower triangular matrix
        t(2:7) = lower_triangular2tensor(t(2:7));
    else
        t = lsqcurvefit(@dtifun, x0, bmat, datafit(:,v), lb, ub, fitoptions);
    end
    
    tensor(:,v) = t(2:7);
    S0_fit(v) = t(1);
end

fprintf('Fit performed in %.2f seconds. \n', toc);

%% rearrange output
sz_fit = sz(2:end);
if length(sz_fit)==1
    sz_fit = [1 sz_fit];
end

dti_fit.S0 = zeros(sz_fit);
dti_fit.S0(sels) = S0_fit;

dti_fit.D = squeeze(zeros([6, sz_fit]));
dti_fit.D(:,sels) = tensor;
dti_fit.D = permute(dti_fit.D, [2:length(size(dti_fit.D)) 1]);
end
