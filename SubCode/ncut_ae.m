% Angular embedding (generalization of normalized cuts).
% Return a set of eigenvectors along with corresponding weights.
function [evecs evals] = ncut_ae(C, O, nvec)
   % get matrix size
   sz = length(C);
   % specify default number of eigenvectors
   if (nargin < 3), nvec = min(sz,8); end
   % specify options
   opts.issym=0;
   opts.isreal = 1;
   opts.disp=2;
   % compute confidence degree matrix
   C1 = sum(C,2);
   C_tot = sum(C1,1);
   % compute normalized total-confidence degree matrix
   D = spdiags(C1.*1./C_tot,0,sz,sz);
   % compute measurement matrix
   MC = (spdiags(1./C1,0,sz,sz)*C);
   [x y v] = find(MC);
   inds = sub2ind([sz sz],x,y);
   MO = sparse(x,y,exp(1i.*full(O(inds))),sz,sz);
   M = MC.*MO;
   % compute error matrix
   I = speye(sz,sz);
   W = ((I - M)')*D*((I - M));
   % compute eigenvectors
   [evecs evals] = eigs(W, D, nvec, 'sm', opts);
   evals = diag(evals);
end
