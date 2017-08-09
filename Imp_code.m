%% Plotting Red green tumor - patch importances
%Xts = lasso_fista(Yts, D, [],opts.lambda, opts);


patches = {};
for i = 8%1:length(tsinds) % Choose the test patient index
    tsfeats = [];
    im = double(imcells{tsinds(i)});
    c = size(im, 3);
    nrow = size(im,1);
    ncol = size(im,2);
    midslice = round(c/2); %Calculate the importance plot from the midslice
    for m = 1 : nrow-psize
        m
        for n = 1: ncol-psize
            patches{m}{n} = im(m:m+psize-1, n:n+psize-1, midslice);
            tsfeats = [tsfeats; patches{m}{n}(:)'];
        end
    end    
end

Yts = tsfeats';
a = DLSI_pred(Yts, D, opts);
pred = []; f1s = [];
a = a-1;
a1=reshape(a,[ncol-psize, nrow-psize]);
a1 = a1';
%subplot(1,2,1); imshow(im(:,:,midslice));
%subplot(1,2,2); imshow(a1);

totmat = zeros(nrow,ncol);
tempmat = zeros(nrow,ncol);
for m = 1 : nrow-psize
    m
    for n = 1: ncol-psize
        totmat(m:m+psize-1, n:n+psize-1) = totmat(m:m+psize-1, n:n+psize-1) + a1(m,n);
        tempmat(m:m+psize-1, n:n+psize-1) = tempmat(m:m+psize-1, n:n+psize-1) + 1;
    end
end
fin_pred = overall_pred(i)
if(fin_pred == 1)
    totmat = totmat./tempmat;
else
    totmat = (tempmat - totmat)./tempmat;
end


temp = im(:,:,midslice);
subplot(1,2,1); imagesc(temp(psize:nrow-psize, psize:ncol-psize)); colorbar
subplot(1,2,2); imagesc(totmat(psize:nrow-psize, psize:ncol-psize)); colorbar
