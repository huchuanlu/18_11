function W = computeW(is_old,imageX,dataW,emag,ephase,edges,sp_ind_middle)
% W = computeW(imageX,dataW,emag,ephase)
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.
[p,q] = size(imageX);

fprintf('Getting indices...\n');


% if is_old==1
%     [w_i,w_j] = cimgnbmap([p,q],dataW.sampleRadius,dataW.sample_rate);
% else
%     [w_i,w_j] = my_cimgnbmap([p,q],dataW.sampleRadius,dataW.sample_rate,dataW.innerRadius);
% end

% sp_ind = unique(edges(:,2));
% % w_i = [0;edges(:,1)];
% w_i = edges(:,1);
% w_j = zeros(length(sp_ind)+1,1);
% w_j(1,1)=0;
% 
% for jj = 1:length(sp_ind)
%     sort_number = edges(1+w_j(jj,1),2);
%     number = sum(edges(:,2)==sort_number);
%     w_j(jj+1,1) = w_j(jj,1)+number;
% end
% % w_j(w_j~=0) = w_j(w_j~=0)+1;
% 
% fprintf('Computing affinities...\n');
% w_i = uint32(w_i);
% w_j = uint32(w_j);
% W = affinityic(emag,ephase,w_i,w_j,sp_ind_middle,max(emag(:)) * dataW.edgeVariance);

w_i = edges(:,1);
w_j = edges(:,2);
w_i = w_i-1;
w_j = w_j-1;
w_i = uint32(w_i);
w_j = uint32(w_j);
W = affinityic_my(emag,ephase,w_i,w_j,sp_ind_middle,max(emag(:)) * dataW.edgeVariance);
%W = affinityic(emag,ephase,w_i,w_j);
%W = my_affinityic(emag,ephase,w_i,w_j);
%W = my_affinityic(emag,ephase,w_i,w_j,max(emag(:)) * dataW.edgeVariance);
fprintf('Done building affinity matrix\n');

%W = W/max(W(:));