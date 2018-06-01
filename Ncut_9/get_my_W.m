function [W,ephase]=get_my_W(im,edge_im,edges,sp_ind_middle)
%function [W1,W2]=get_my_W(im,edge_im)


edge_im=double(edge_im);


%blurring 
G = fspecial('gaussian',[15 15],3);      
edge_im = imfilter(edge_im,G,'same');

[W,imageEdges,ephase] = ICgraph(rgb2gray(im),edge_im,edges,sp_ind_middle);
%[W1,W2,imageEdges] = my_ICgraph(rgb2gray(im),edge_im);

end

