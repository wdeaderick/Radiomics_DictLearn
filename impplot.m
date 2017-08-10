%% Plotting Red green tumor - patch importances
 %Xts = lasso_fista(Yts, D, [],opts.lambda, opts);
load wspace.mat
load dict_idh1.mat
psize = 16;
totalpat = 140;
patches = {};
 for i = 94  % Choose the test patient index
     tsfeats = [];
     im = double(imcells_seq4{i});
     c = size(im, 3);
     nrow = size(im,1);
     ncol = size(im,2);
     midslice = round(c/2); %Calculate the importance plot from the midslice
     for m = 1 : nrow-psize
         for n = 1: ncol-psize
             patches{m}{n} = im(m:m+psize-1, n:n+psize-1, midslice);
             tsfeats = [tsfeats; patches{m}{n}(:)'];
         end
     end
 end
 
Yts = tsfeats';
k = 300;   %Number of dictionary atoms
C              = 2;
D_range        = k*(0:C);
opts.k = k;
opts.lambda    = 0.1;
opts.eta       = 1;
opts.D_range   = D_range;
opts.show_cost = 0;
opts.show      = 0;
opts.verbose    = false;
opts.max_iter  = 10;
a = DLSI_pred(Yts, dict_1p19q, opts);
pred = []; f1s = [];
a = a-1;
 
for i = 1:size(Yts,2)/totalpat
    lookvec = [(i-1)*totalpat+1:i*totalpat];
    [temp] = mode(a(lookvec));
    f1s = [f1s; sum(a(lookvec) == 1 )]; %Confidence of prediction
    pred = [pred; temp];
end
pred
clabels(96)
pred = mode(pred);
a1=reshape(a,[ncol-psize, nrow-psize]);
a1 = a1';
%subplot(1,2,1); imshow(im(:,:,midslice));
%subplot(1,2,2); imshow(a1);
 
totmat = zeros(nrow,ncol);
tempmat = zeros(nrow,ncol);
for m = 1 : nrow-psize
    for n = 1: ncol-psize
        totmat(m:m+psize-1, n:n+psize-1) = totmat(m:m+psize-1, n:n+psize-1) + a1(m,n);
        tempmat(m:m+psize-1, n:n+psize-1) = tempmat(m:m+psize-1, n:n+psize-1) + 1;
    end
end
if(pred == 1)
    totmat = totmat./tempmat;
else
    totmat = (tempmat - totmat)./tempmat;
end
 
 
temp = im(:,:,midslice);
subplot(1,2,1); imagesc(temp(psize:nrow-psize, psize:ncol-psize)); colorbar; title('Original Image')
subplot(1,2,2); imagesc(totmat(psize:nrow-psize, psize:ncol-psize)); colorbar; title('PI Plot')
