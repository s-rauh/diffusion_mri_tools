%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates MD, FA and Eigenvalues of a diffusion tensor.
% The tensor elements should be given in the following form: 
% dimensions data (= tensor): (x,y,z,tensor-elements)
% tensor-elements(1) = Dxx
% tensor-elements(2) = Dyy
% tensor-elements(3) = Dzz
% tensor-elements(4) = Dxy
% tensor-elements(5) = Dxz
% tensor-elements(6) = Dyz
%
% Returns a structure dti_params with the following fields:
%   - MD: Mean diffusivity
%   - FA: Fractional anisotropy
%   - eigenval: The 3 eigenvalues of the diffusion tensor
%   - eigenvec: The corresponding 3 eigenvectors of the diffusion tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dti_params = calc_dti_parameters(tensor)

%reshape tensor
if size(tensor, ndims(tensor)) ~= 6
    error('Tensor elements expected in last dimension')
end

tensor = permute(tensor, [length(size(tensor)), 1:length(size(tensor))-1]);
tensor_sz = size(tensor);
tensor_dim = ndims(tensor);

tensor = reshape(tensor, 6, []);

MD = zeros(size(tensor,2),1);
FA = zeros(size(tensor,2),1);
eval = zeros([size(tensor,2),3], 'single');
evec = zeros([size(tensor,2),3,3], 'single');

fprintf('Calculate pyqmri DTI parameters for each voxel... \n');
tic;
for x = 1:size(tensor,2)
    temp = squeeze(tensor(:,x));
    %put tensor elements together
    dti = [temp(1) temp(4) temp(5); temp(4) temp(2) temp(6); temp(5) temp(6) temp(3)];
    %calculate eigenvectors, eigenvalues, FA and MD for each voxel
    [Eigenvectors, D] = eig(dti);
    EigenValues = diag(D);
    [~, index] = sort(EigenValues, 'descend');
    EigenValues = EigenValues(index); Eigenvectors = Eigenvectors(:,index);
    
    MDv = (EigenValues(1)+EigenValues(2)+EigenValues(3))/3;
    FAv = sqrt(1.5)*(sqrt((EigenValues(1)- MDv).^2 + (EigenValues(2) - MDv).^2 + (EigenValues(3) - MDv).^2) ...
        ./ sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2));
    
    % store MD, FA and Eigenvalues
    MD(x) = MDv;
    FA(x) = FAv;
    eval(x,:) = EigenValues;
    evec(x,:,:) = Eigenvectors;
end

%reshape to matrix size
dti_params.MD = reshape(MD, tensor_sz(2:tensor_dim));
dti_params.FA = reshape(FA, tensor_sz(2:tensor_dim));
dti_params.eigenval = reshape(eval, [tensor_sz(2:tensor_dim), 3]);
dti_params.eigenvec = reshape(evec, [tensor_sz(2:tensor_dim), 3, 3]);

fprintf('Tensor calculation finished in %.2f seconds. \n', toc);
end
