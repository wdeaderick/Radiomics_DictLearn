function [D, pred, f1s] = runCOPAR(Y, Yrange, Yts, totalpat)
%Change COPAR input parameters here
fprintf('Running COPAR... \n')
k = 20;
k0 = 10;
opts.k = k;
opts.k0 = k0;
C                = 2;
D_range_ext = [k*(0:C), k*C + k0];
opts.lambda      = 0.001;
opts.eta         = 0.01;    
opts.show        = false;
opts.max_iter    = 10;        
opts.verbose      = false;

tic
[D, X, rt] = COPAR(Y, Yrange, opts);
toc
opts.classify_mode = 'LC'; %'LC' or 'GC'
opts.gamma = 0.01; %[0.0001, 0.001, 0.005, 0.01]
a = COPAR_pred(Yts, D, D_range_ext, opts);

pred = []; f1s = [];
a = a-1;
for i = 1:size(Yts,2)/totalpat
    lookvec = [(i-1)*totalpat+1:i*totalpat];
    [temp] = mode(a(lookvec));
    f1s = [f1s; sum(a(lookvec) == 1 )]; %Confidence of prediction
    pred = [pred; temp];
end
