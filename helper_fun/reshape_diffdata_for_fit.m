%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to reshape diffusion data to a bval x n array for data fitting.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%   -   data
%   -   b-values
%
% Output:
%   -   diffdata
%       Reshaped diffusion data
%       The first dimension is the b-value and the second dimension 
%       contains the data. 
%   -   data_sz
%       The original array size (data_sz) is returned to allow reshaping 
%       of the data. 
%   -   Optional:
%       The b-value dimension can be returned (varargout{1}).
%       Permutation vector can be returned (varargout{2}).         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [diffdata, data_sz, varargout] = reshape_diffdata_for_fit(data, bval)

try
    bdim = find(size(data)==length(bval));
catch
    error('Datasize and number of b-values are not matching!')
end
diffdata = permute(data, [bdim, 1:bdim-1 bdim+1:length(size(data))]);
data_sz = size(diffdata);

diffdata = reshape(diffdata, length(bval), []);
diffdata = double(diffdata);

if nargout > 2
    varargout{1} = bdim;
end

if nargout > 3
    %permutation vector
    varargout{2} = [bdim, 1:bdim-1 bdim+1:length(size(data))];
end