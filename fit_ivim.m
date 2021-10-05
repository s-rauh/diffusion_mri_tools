%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to perform an IVIM fit to data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input values: 
%   - Diffusion data:   Array containing the data. Can be single voxel, a
%                       slice or volumetric data. The last dimension needs
%                       to match the number of b-values. 
%
%   - b-values:         Vector containing the b-values
%
% Options:
%   - initialguess:     S0, D, f, Ds
%                       If not provided, default guess will be used.
%
%   - mask:             Default 1
%                       Background is masked by using a threshold cut-off
%                       value.
%
%   - normalize:        Default 1
%                       Normalize data to the min(bval)-signal
%
%   - fit_method:       Default 'two_step'
%                       'free', 'two_step', 'segmented'
%                       'free': The full data with all b-values is fitted 
%                       to the IVIM equation. 
%                       'two_step': A two-step fit will be performed,                   
%                       estimating D and f from high b-values (b >= bcut) 
%                       only and performing a bi-exponential fit with fixed 
%                       D to obtain f and D*.
%                       'segmented': A segmented IVIM fit is performed, 
%                       using b = 0 and b >= bcut to estimate D and f. 
%                       D* is not fitted with this method. 
%
%   - bcut:             Default 200
%                       Cut-off b-value for the two_step and segmented fit.
%
%   - seg_data:         Default 0
%                       Select only data with bval >= bcut and b = 0 for
%                       fitting. Only used for segmented IVIM fit
%                       (fit_method = 'segmented')
%
% Output:
% Structure ivim_fit containing the following fields:
%   - S0 signal (fitted). In case of normalization S0 will be close to 1
%     for all voxels.
%   - diffusion coefficient D
%   - perfusion fraction f
%   - pseudo-diffusion coefficient Ds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ivim_fit = fit_ivim(data, bval, options)

arguments
    data
    bval
    options.initialguess (4,1)
    options.mask {mustBeNumericOrLogical} = 1
    options.normalize {mustBeNumericOrLogical} = 1    
    options.fit_method {mustBeMember(options.fit_method, {'free', 'two_step', 'segmented'})} = 'two_step'
    options.bcut {mustBeNumeric} = 200
    options.seg_data {mustBeNumericOrLogical} = 0

end

% set initial guess for fit
if isfield(options, 'initialguess')
    x0 = options.initialguess;
else
    x0 = [1, 1.7, 0.1, 10];
end

%% rearrange and reshape data to a bval x n array
[data, sz] =  reshape_diffdata_for_fit(data, bval);

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
    x0(1) = 1;  %initial guess S0
else
    x0(1) = mean(datafit(bval==min(bval),:), 'all');
end

%% scale b-value
%scale b-value
if size(bval,1)==1
    bval = bval.';
end
bval = bval_scaling(bval);
options.bcut = bval_scaling(options.bcut);
%% set fit options
%boundaries
lb = [0,   0, 0, 5];
ub = [inf, 5, 1, 300];

if options.normalize
    ub(1) = 10;
end

fitoptions = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt', ...
    'FunctionTolerance', 1e-10, 'MaxIterations', 1000, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
    'Display', 'off');

%initialize fit parameters
S0_fit = zeros(1, size(datafit,2));
D_fit = zeros(1, size(datafit,2));
f_fit = zeros(1, size(datafit,2));
Ds_fit = zeros(1, size(datafit,2));
%% fit data voxelwise
fprintf('Perform voxel-wise IVIM fit...\n')
tic;
for v = 1:size(datafit, 2)
    switch options.fit_method
        %%%%%%%%%%%%%%%%%%%
        % free fit
        %%%%%%%%%%%%%%%%%%
        case 'free'
            x = lsqcurvefit(@ivimfun, x0, bval, datafit(:,v), lb, ub, fitoptions);
        
        %%%%%%%%%%%%%%%%%%%
        % two_step fit
        %%%%%%%%%%%%%%%%%%
        case 'two_step'
            if v == 1
                %adjust boundaries and initial guess
                lb_twostep = lb([1 3 4]); ub_twostep = ub([1 3 4]); x0_twostep = x0([1 3 4]);
            end
            
            %First, perform mono-exponential fit to b >= bcut to get an estimate for D
            [param(1:2)] = lsqcurvefit(@(x, xdat) x(1)*exp(-xdat*x(2)), x0([1,2]), bval(bval>=options.bcut), datafit(bval>=options.bcut,v), ...
                lb([1,2]), ub([1,2]), fitoptions);
            
            x(2) = param(2);    %D
            
            %Estimate f and use it as initial guess for IVIM fit
            S0_tmp = mean(datafit(bval==min(bval),v));
            f_guess = 1 - (param(1) / S0_tmp);
            if f_guess < 0 || isnan(f_guess)
                f_guess = 0;
            end
            x0_twostep(2) = f_guess;
            
            %perform fit. D is fixed for the bi-exponential IVIM fit.
            input.bval = bval;
            input.D_fix = x(2);
            x([1,3:4]) = lsqcurvefit(@ivimfun, x0_twostep, input, datafit(:,v), lb_twostep, ub_twostep, fitoptions);
        
        %%%%%%%%%%%%%%%%%%%
        % segmented fit
        %%%%%%%%%%%%%%%%%%
        case 'segmented'
            %optional: select desired b-value range (if full range was
            %acquired)
            if options.seg_data && v == 1
               datafit = [datafit(bval==0,:); datafit(bval>=options.bcut,:)];
               bval = [bval(bval==0); bval(bval>=options.bcut)];
            end
            
            %perform fit
            x = lsqcurvefit(@ivimfun_seg, x0(1:3), bval, datafit(:,v), lb(1:3), ub(1:3), fitoptions);
            
    end
    
    S0_fit(v) = x(1);
    D_fit(v) = x(2);
    f_fit(v) = x(3);
    if ~isequal(options.fit_method, 'segmented')
        Ds_fit(v) = x(4);
    end
end

fprintf('Fit performed in %.2f seconds. \n', toc);

%% rearrange output
sz_fit = sz(2:end);
if length(sz_fit)==1
    sz_fit = [1 sz_fit];
end
ivim_fit.S0          = zeros(sz_fit);
ivim_fit.S0(sels)    = S0_fit;

ivim_fit.D           = zeros(sz_fit);
ivim_fit.D(sels)     = D_fit;

ivim_fit.f           = zeros(sz_fit);
ivim_fit.f(sels)     = f_fit;

if ~isequal(options.fit_method, 'segmented')
    ivim_fit.Ds          = zeros(sz_fit);
    ivim_fit.Ds(sels)    = Ds_fit;
end

end
