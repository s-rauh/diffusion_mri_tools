%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to perform a combined IVIM-DTI fit to data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input values: Diffusion data, b-values, diffusion gradient directions
% Options:
%   - mask:             Default 1
%                       Background is masked by using a threshold cut-off
%                       value based on the minimum b-value data. 
%
%   - normalize:        Default 1
%                       Normalize data to the min (bval)-signal
%
%   - initialguess:     Initial guess for f and Ds
%                       If not provided, default guess will be used. 
%                       Only f and Ds can be set here, the initial guess
%                       for S0 and the diffusion tensor is kept constant. 
%
%   - dti_constrained   Default 1
%                       Constrained fit for the diffusion tensor elements.
%                       The lower triangular matrix of a Cholesky
%                       decomposition of the diffusion tensor is fitted.
%                       After the fit, the diffusion tensor can be
%                       calculated. The fit boundaries differ from the
%                       unconstrained fit. 
%
%   - fit_method        'free', 'two_step', 'IVIM_correct'
%                       'free': The data is fitted to the full IVIM-DTI
%                       equation. 
%                       'two_step': First, the diffusion tensor and f are
%                       estimated using the high b-values only (b >= bcut).
%                       Then, Dtensor is kept constant and a bi-exponential
%                       fit is performed to estimate f and Ds. The
%                       previously calculated value for f is used as
%                       initial guess. 
%                       'IVIM_correct': First an IVIM fit is performed to
%                       obtain f and Ds. The data is then IVIM-corrected by
%                       substracting the IVIM component. A pure DTI fit is
%                       performed on the remaining data. 
%
%   - bcut:             Default 200
%                       Cut-off b-value for the two_step fit. 
%
% Output:
% Structure ivimdti_fit containing the following fields:
%   - S0 signal (fitted). In case of normalization S0 will be close to 1
%     for all voxels. 
%   - tensor: Diffusion tensor containing the 6 elements in the order
%             Dxx, Dyy, Dzz, Dxy, Dxz, Dyz
%   - perfusion fraction f
%   - pseudo-diffusion coefficient Ds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ivimdti_fit = fit_ivimdti(data, bval, diffdir, options)

arguments
    data 
    bval
    diffdir
    options.mask {mustBeNumericOrLogical} = 1
    options.normalize {mustBeNumericOrLogical} = 1
    options.initialguess (2,1)
    options.dti_constrained = 1;
    options.fit_method {mustBeMember(options.fit_method, {'free', 'two_step', 'IVIM_correct'})} = 'free'
    options.bcut {mustBeNumeric} = 200
end

%rearrange and reshape data to a bval x n array
[data, sz] =  reshape_diffdata_for_fit(data, bval);

%% optional: mask background voxels
%define datafit
if options.mask
    [datafit, sels] = mask_diffdata(data, bval);
else
    sels = true(1, size(data,2));
    datafit = data;
end

%% optional: normalize data
if options.normalize
    datafit = norm_diffdata(datafit, bval);
    %set initial guess
    x0 = [1 1 1 1 0 0 0 0.1 10];
else
    x0 = [mean(datafit(bval==min(bval),:), 'all') 1 1 1 0 0 0 0.1 10];
end

%% scale b-value
if size(bval,1)==1
    bval = bval.';
end
if max(bval) > 100
    bval = bval./1000;
end
if options.bcut > 1
    options.bcut = options.bcut/1000;
end

%% set fit options
diffparams.bval = bval;
diffparams.diffdir = diffdir;

if isfield(options, 'initialguess')
    %initial guess for f and D*
    x0(8:9) = options.initialguess;
end

if options.dti_constrained
    %Fit Cholesky lower triangular matrix coefficients
    lb = [0   -9 -9 -9 -9 -9 -9 0   5];
    ub = [inf  9  9  9  9  9  9 1 300];
else
    lb = [0   0 0 0 0 0 0 0   5];
    ub = [inf 3 3 3 3 3 3 1 300];
end
if options.normalize
    ub(1) = 10;
end

fitoptions = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt', ...
    'FunctionTolerance', 1e-10, 'maxIterations', 1000, 'OptimalityTolerance', 1e-10, ...
    'Display', 'off');

%initialize fit parameters
S0_fit = zeros(1, size(datafit, 2));
tensor_fit = zeros(6, size(datafit, 2));
f_fit = zeros(1, size(datafit, 2));
Ds_fit = zeros(1, size(datafit, 2));

%% perform fit
fprintf('Perform voxel-wise IVIM-DTI fit...\n');
tic;
%Select fit method
switch options.fit_method
    
    %%%%%%%%%%%%%%%%%%%
    % free fit
    %%%%%%%%%%%%%%%%%%
    case 'free'
        %loop over voxels
        for v = 1:size(datafit,2)
            if options.dti_constrained
                x = lsqcurvefit(@ivimdtifun_constr, x0, diffparams, datafit(:,v), lb, ub, fitoptions);
                x(2:7) = lower_triangular2tensor(x(2:7));
            else
                x = lsqcurvefit(@ivimdtifun, x0, diffparams, datafit(:,v), lb, ub, fitoptions);
            end
            S0_fit(:,v) = x(1);
            tensor_fit(:,v) = x(2:7);
            f_fit(:,v) = x(8);
            Ds_fit(:,v) = x(9);
        end
        
    %%%%%%%%%%%%%%%%%%%
    % two-step fit
    %%%%%%%%%%%%%%%%%%
    case 'two_step'
        %loop over voxels
        for v = 1:size(datafit,2)
            if v == 1
                lb_seg = lb([1 8 9]); ub_seg = ub([1 8 9]); x0_seg = x0([1 8 9]);
            end
            %first fit tensor to b >= bcut
            if options.dti_constrained
                [param(1:7)] = lsqcurvefit(@dtifun_constr, x0(1:7), calc_bmat(bval(bval>=options.bcut), diffdir(bval>=options.bcut,:)), ...
                    datafit(bval>=options.bcut, v), lb(1:7), ub(1:7), fitoptions);
                x(2:7) = lower_triangular2tensor(param(2:7));
            else
                [param(1:7)] = lsqcurvefit(@dtifun, x0(1:7), calc_bmat(bval(bval>=options.bcut), diffdir(bval>=options.bcut,:)), ...
                    datafit(bval>=options.bcut,v), lb(1:7), ub(1:7), fitoptions);
                x(2:7) = param(2:7);
            end
            %estimate f
            S0_tmp = mean(datafit(bval==min(bval), v));
            f_guess = 1 - (S0_tmp / param(1));
            if f_guess < 0
                f_guess = 0;
            end
            x0_seg(2) = f_guess;
            
            %perform fit with fixed diffusion tensor
            diffparams.tensor = x(2:7);
            x([1,8:9]) = lsqcurvefit(@ivimdtifun, x0_seg, diffparams, datafit(:,v), lb_seg, ub_seg, fitoptions);
            
            
            S0_fit(:,v) = x(1);
            tensor_fit(:,v) = x(2:7);
            f_fit(:,v) = x(8);
            Ds_fit(:,v) = x(9);
        end
        
    %%%%%%%%%%%%%%%%%%%
    % IVIM-correct DTI data
    %%%%%%%%%%%%%%%%%%    
    case 'IVIM_correct'      
        %fit IVIM to the signal (ivimfun)
        ivim_pars = fit_ivim(datafit, bval, 'mask', 0, 'normalize', 0, 'two_step', 1);
        %IVIM-correct DTI signal by substracting IVIM component
        datafit_cor = datafit - ivim_pars.S0.*ivim_pars.f.*exp(-bval.*ivim_pars.Ds);
        %fit DTI to remaining DTI signal (dtifun)
        dti_pars = fit_dti(datafit_cor, bval, diffdir, ...
            'constrained', options.dti_constrained, 'mask', 0, 'normalize', 0);      
        
        S0_fit = ivim_pars.S0;
        f_fit = ivim_pars.f;
        Ds_fit = ivim_pars.Ds;
        tensor_fit = dti_pars.D.';
end
fprintf('Fit performed in %.2f seconds. \n', toc);

%% rearrange output
sz_fit = sz(2:end);
if length(sz_fit)==1
    sz_fit = [1 sz_fit];
end

ivimdti_fit.S0       = zeros(sz_fit);
ivimdti_fit.S0(sels) = S0_fit;

ivimdti_fit.tensor           = squeeze(zeros([6, sz_fit]));
ivimdti_fit.tensor(:,sels)   = tensor_fit;
ivimdti_fit.tensor           = permute(ivimdti_fit.tensor, ...
    [2:length(size(ivimdti_fit.tensor)) 1]);

ivimdti_fit.f        = zeros(sz_fit);
ivimdti_fit.f(sels)  = f_fit;

ivimdti_fit.Ds       = zeros(sz_fit);
ivimdti_fit.Ds(sels) = Ds_fit;
