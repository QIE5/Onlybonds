function stack_filt = svdClutterFilter(frame,nCut,nNoiseCut)

%Flatten each frame, stack them into a Casorati matrix X (pixels × time)
%perform SVD,remove low-rank (tissue clutter) or high-rank (noise) components
%and then reconstruct the frames
frame = im2double(frame); %linear calculation so use double
[H,W,T]=size(frame);
S=H*W;
X=reshape(frame,S,T); %3D sequence to 2D Casorati matrix X
%if S>>T then use temporal covariance
C = (X'*X); %dimension is T*T
% Eigen-decomposition: C=V*D*V' and S*S=VΣ^2V^*
[V, D] = eig(C, 'vector'); %eig(c)to return eigenvector and eigenvalue, use vector to return vector instead of diagonal matrix
%sorting, because we need to eliminate the first few items with largerst energy 
[Dsorted, idx] = sort(real(D), 'descend');%idx=index map
V = V(:, idx); %column indexing to make V line up with D


sigma = sqrt(max(Dsorted, 0));%or use sigma =sqrt(Dsorted)
epsSigma = 1e-12; %make a lower limit so no NaN or Inf occur
invSigma = 1./max(sigma, epsSigma); %1/sigma use./ because Element wise division
U = X * (V .* invSigma.'); %most important part:Compute U implicitly: U = S * V * inv(Sigma)


%use keep to remove first few and last singular value
keep = true(T,1);
keep(1:nCut) = false;% remove tissue (largest singular values)
%remove noise tail (smallest singular values)
if nNoiseCut > 0
    keep(end-nNoiseCut+1:end) = false;
end


Uf = U(:, keep); %use keep function to eliminate unwanted data
Vf = V(:, keep);
sigmaf = sigma(keep);
Sf=Uf * diag(sigmaf)*Vf'; %or use Uf * (sigmaf .* (Vf'))


stack_filt = reshape(Sf, [H,W,T]);%reshape back
end

