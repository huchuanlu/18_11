function Saliency_AMC_AE
%%  setup parameters and paths
clear;     tic;
addpath( './SubCode/' );
addpath( './Ncut_9/' );
theta=10; % control the edge weight 
tao = 0.3; % threshhold of AE
KK=[250,150,350]; % number of superpixel
imgRoot = './image/';
FCNfeatureRoot = './FCN-feature/';
marginRoot='./boundary-map/';
maskNames=dir([marginRoot '*' '.png']);
SalmapRoot = './saliency-map/';
mkdir( SalmapRoot );
ImgNames=dir([imgRoot '*' 'jpg']);  
ImgNum=length(ImgNames);

%% 
for i=1:ImgNum
    i    
    imname = [ imgRoot ImgNames(i).name ];
    Img = double( imread( imname ) );
    [ height,width ] = size(Img(:,:,1));
    PixNum = height*width;    
    ImgVecR = reshape( Img(:,:,1)', PixNum, 1);
    ImgVecG = reshape( Img(:,:,2)', PixNum, 1);
    ImgVecB = reshape( Img(:,:,3)', PixNum, 1);
    
    m = 20;  % m is the compactness parameter, k is the super-pixel number in SLIC algorithm    
    ImgAttr=[ height ,width, KK(1), m, PixNum ];
    % obtain superpixel from SLIC algorithm: LabelLine is the super-pixel label vector of the image,
    % Sup1, Sup2, Sup3 are the mean L a b colour value of each superpixel,
    % spnum is the number of the super-pixel.
    [ LabelLine, Sup1, Sup2, Sup3, spnum ] = SLIC( ImgVecR, ImgVecG, ImgVecB, ImgAttr );
    Label=reshape(LabelLine,width,height);
    Label = Label';     % the superpixle label
    Label = Label+1;
    [ ConPix, ConPixDouble ] = find_connect_superpixel_DoubleIn_Opposite( Label-1, spnum, height ,width );
    
    % get edges
    edges = [];
    for j=1:spnum
        for z=j+1:spnum
            if ConPixDouble(j,z)>0
                edges = [edges;[j,z]];          
            end
        end
    end
      
    % extract FCN feature
    inds=cell(spnum,1);
    for ii=1:spnum
        inds{ii}=find(Label==ii);
    end   
    [meanVgg1,meanVgg2] = ExtractFCNfeature(FCNfeatureRoot,ImgNames(i).name(1:end-4),inds,height,width);

    % compute affinity matrix
    weights1 = makeweights(edges,meanVgg1,theta);    
    weights2 = makeweights(edges,meanVgg2,theta);  
    W1 = adjacency(edges,weights1,spnum); 
    W2 = adjacency(edges,weights2,spnum);
    dd = sum(W1); D1 = sparse(1:spnum,1:spnum,dd); clear dd;
    L1 =D1-W1; 
    dd = sum(W2); D2 = sparse(1:spnum,1:spnum,dd); clear dd;
    L2 =D2-W2; 
    WL = ComputeAffinityMatrix(L1,L2,spnum);      % the learnt affinity matrix 
  
    WconFirst = diag(diag(WL))\WL;
    Discard = sum(WconFirst,2);
    DiscardPos = find( Discard < 1.1 );    % to discard the outlier
    LenDiscardPos = length(DiscardPos);

    EdgSup = Find_Edge_Superpixels( Label-1, spnum,  height, width , WconFirst, ConPix );    
    for j=1:LenDiscardPos
        EdgSup( DiscardPos(j) ) = 2;
    end
    NumIn = spnum - length( find( EdgSup == 2 ) );    
    Edge_ind = (EdgSup==1);
    EdgWcon = WconFirst(:,Edge_ind);

%% Absorbing Markov Chain     
    alph = 1;  
    Edge_ind = (EdgSup<2);
    BaseEdg = sum( EdgWcon(Edge_ind, : ),2 );
    sumD = Discard(Edge_ind);      
    D = diag( BaseEdg + sumD );
    Wmid = WconFirst(Edge_ind,Edge_ind);
    Wmid = D \ Wmid;
    I = eye( NumIn );
    N = ( I - alph* Wmid  );
    y = ones( NumIn, 1 );
    Sal = N \ y;
    Sal = normalize(Sal);
    SalAll = zeros(spnum,1);           
    SalAll(Edge_ind) = Sal;
    for j=1:LenDiscardPos
        for z=1:spnum
            if ConPix( DiscardPos(j), z ) > 0
                if SalAll(z) >.3
                    SalAll( DiscardPos(j) ) = 1 ;
                    break;
                end
            end
        end
    end
    SalLine = sup2pixel( PixNum, LabelLine, SalAll );  % to convey the saliency value from superpixel to pixel 
    Salpix = reshape( SalLine, width, height );
    Salpix = Salpix';
    imwrite( Salpix, [ SalmapRoot, ImgNames(i).name(1:end-4),sprintf('_WL') '.png' ] );

%% Angular Embedding
     imname = [ SalmapRoot, ImgNames(i).name(1:end-4),sprintf('_WL') '.png' ];
     sal_map = double( imread( imname ) );   
     sal_map1 = zeros(spnum,1);
     for ii=1:spnum
          sal_map1(ii)=sum(sal_map(inds{ii}))./length(inds{ii});    
     end
     sal_map1 = normalization(sal_map1, 0); 
     thresh_H = 0.95;
     seed_fore = sal_map1>thresh_H;
     thresh_L = 0.1;   
     seed_back= sal_map1<thresh_L;
     seed_middle = (~seed_fore).*(~seed_back)>0.5;
     my_edges = GetEdges(seed_middle,seed_middle) ;
     if isempty(my_edges)
          outname1=[SalmapRoot ImgNames(i).name(1:end-4) sprintf('_WL_AE')  '.png'];    
          imwrite(Salpix,outname1);
          continue;
     end

     % boundary map
     maskName=[marginRoot maskNames(i).name];
     mask = imread(maskName);
     mask_new = im2double(mask);
     sp_ind_middle = zeros(spnum,1);
     for ii=1:spnum
         sp_ind_middle(ii)=median(inds{ii});    
     end

     sp_ind_middle = uint32(sp_ind_middle);
     [W_edges,~]=get_my_W(Img,mask_new,my_edges,sp_ind_middle);
     W_edges=full(W_edges);
     C_AE1 = full(WL);
     C_AE = C_AE1(seed_middle,seed_middle);     
     O = sal_map1(my_edges(:,1))-sal_map1(my_edges(:,2)); 
     O = O.*(1~=W_edges);
     O = O.*(abs(O)>tao);
     O = sign(O)*pi/8;
     O_sum = sum(abs(O));
     if O_sum==0
         outname1=[SalmapRoot ImgNames(i).name(1:end-4) sprintf('_WL_AE')  '.png'];   
         imwrite(Salpix,outname1);
         continue;
      end
      O=O/O_sum*(pi/2); 
      O_full = adjacency_skew_symmetric(my_edges,O,spnum);
      O_full = full(O_full);
      O_AE = O_full(seed_middle,seed_middle);

     [evecs,~] = ncut_ae(C_AE, O_AE, 1);
     ordering = angle(evecs);
     if length(unique(ordering))==1
         outname1=[SalmapRoot ImgNames(i).name(1:end-4) sprintf('_WL_AE')  '.png']; 
         imwrite(Salpix,outname1);
         continue;
     end
     
     % --------------------- 对中间灰度的显著性归一化 ------------------------  
     ordering = normalization(ordering, 0); 
     sal_map = zeros(spnum,1);
     sal_map(seed_middle) = ordering;
     sal_map(seed_fore) = 1;        %将大于thresh_H的部分强制设置为1，保证准确性
     % output salmap 
     tmapstage2=zeros(height,width);
     for ii=1:spnum
         tmapstage2(inds{ii})=sal_map(ii);    
     end
        
     tmapstage2=normalization(tmapstage2, 0);
     outname1=[SalmapRoot ImgNames(i).name(1:end-4) sprintf('_WL_AE')  '.png'];   
     imwrite(tmapstage2,outname1);

end
toc;

