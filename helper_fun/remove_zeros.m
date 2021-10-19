function [data, diffparams] = remove_zeros(data, diffparams)

if isstruct(diffparams)
    %IVIM-DTI
    diffparams.bval(data==0) = [];
    diffparams.diffdir(data==0,:) = [];
elseif isvector(diffparams)
    %IVIM, diffparams = bval
    diffparams(data==0) = [];
else
    %DTI, diffparams = b-matrix
    diffparams(data==0,:) = [];
end
data(data==0) = [];

