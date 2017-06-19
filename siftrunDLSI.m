function [D, pred, f1s, Xts] = siftrunDLSI(Y, Yrange, Yts, tslengths)
%Used when employing the SIFT desciptor
%Change DLSI input parameters here
fprintf('Running DLSI... \n')
k = 20;   %Number of dictionary atoms
C              = 2;
D_range        = k*(0:C);
opts.lambda    = 0.001;
opts.eta       = 0.01;
opts.D_range   = D_range;
opts.show_cost = 0;
opts.show      = 0;
opts.verbose    = false;
opts.max_iter  = 10; 
tic
[D, X, rt] = DLSI(Y, Yrange, opts);
toc
Xts = lasso_fista(Yts, D, [],opts.lambda, opts);
a = DLSI_pred(Yts, D, opts);
pred = []; f1s = [];
a = a-1;
tsbookmark = 0;
for i = 1:length(tslengths)
    lookvec = tsbookmark + 1:tsbookmark + tslengths(i);
    [temp] = mode(a(lookvec));
    f1s = [f1s; sum(a(lookvec) == 1 )]; %Confidence of prediction
    pred = [pred; temp];
    tsbookmark = tsbookmark + tslengths(i);
end
 
%{
% Use to code chunk to view Dictionary atoms
m = size(D,2)/10;
for i = 1:size(D,2)
    figure(1);
    subplot(10,m,i);
    patch = reshape(D(:,i), [psize, psize]);
    imtemp = (patch - min(patch(:)))/(max(patch(:)) - min(patch(:)));
    imshow(imtemp, 'InitialMagnification', 2000);
end
%}
