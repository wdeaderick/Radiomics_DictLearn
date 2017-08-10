function [mean_aucs] = pixelensemble()
%The method argument should be passed 1 for DLSI, 2 for COPAR, 3 for FDDL.
%Returns a 1x3 vectors "mean_aucs" containing the mean AUCs for IDH Status, Grade, and Codeletion, respectively.
load wspace.mat
mean_aucs = {};%Preallocate return cell
%allinds is the indices of the 108 patients who have flair, t1, and t2
%MRIs available.
 
seq1masks = {allimagesmasks{1}{allinds}};
seq2masks = {allimagesmasks{2}{allinds}};
seq4masks = {allimagesmasks{4}{allinds}};
 
labels = labels(allinds);
glabels = glabels(allinds);
clabels = clabels(allinds);
hlabels = hlabels(allinds);
slabels = slabels(allinds);
allpws = allpws(allinds,:);
 
 
%% All patients training
labelsall = [labels, glabels, clabels];
for labs = 1
    finlabels = labelsall(:,labs);
    totalpat = 140; %number of patches per patient
    allfeats = {};
    allfeats{1} = {}; %contains all the image patches 
    allfeats{2} = {};
    allfeats{3} = {};
    psize = 16; %patch size: [psize x psize]
    for p = 1:length(finlabels)
        images = {};
        images{1} = double(seq1imcells_ensemble{p});
        images{2} = double(seq2imcells_ensemble{p});
        images{3} = double(seq4imcells_ensemble{p});
        nrow = [];
        ncol = [];
        nrow(1) = size(images{1},1);
        ncol(1) = size(images{1},2);
        c1 = size(images{1}, 3);
        arr = {};
        arr{1} = [];
        arr{2} = [];
        arr{3} = [];
        nrow(2) = size(images{2},1);
        ncol(2) = size(images{2},2);
        c2 = size(images{2}, 3);
        nrow(3) = size(images{3},1);
        ncol(3) = size(images{3},2);
        c4 = size(images{3}, 3);
        for sel = 1:c1
            imtemp = seq1masks{p}(:,:,sel);
            if(sum(imtemp(:)) > 0.50*nrow(1)*ncol(1)) %select slices whose tumor part is more than 30%
                arr{1} = [arr{1} sel];
            end
        end
        if (isempty(arr{1}))
            for sel = 1:c1
                imtemp = seq1masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.30*nrow(1)*ncol(1)) %select slices whose tumor part is more than 30%
                    arr{1} = [arr{1} sel];
                end
            end
        end
        if (isempty(arr{1}))
            for sel = 1:c1
                imtemp = seq1masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.20*nrow(1)*ncol(1)) %select slices whose tumor part is more than 30%
                    arr{1} = [arr{1} sel];
                end
            end
        end
        if (isempty(arr{1}))
            for sel = 1:c1
                imtemp = seq1masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.10*nrow(1)*ncol(1)) %select slices whose tumor part is more than 30%
                    arr{1} = [arr{1} sel];
                end
            end
        end
        if (isempty(arr{1}))
            for sel = 1:c1
                imtemp = seq1masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.04*nrow(1)*ncol(1)) %select slices whose tumor part is more than 30%
                    arr{1} = [arr{1} sel];
                end
            end
        end 
        if (isempty(arr{1}))
            for sel = 1:c1
                imtemp = seq1masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0*nrow(1)*ncol(1)) %select slices whose tumor part is more than 30%
                    arr{1} = [arr{1} sel];
                end
            end
        end
        for sel = 1:c2
            imtemp = seq2masks{p}(:,:,sel);
            if(sum(imtemp(:)) > 0.50*nrow(2)*ncol(2)) %select slices whose tumor part is more than 30%
                arr{2} = [arr{2} sel];
            end
        end
        if (isempty(arr{2}))
            for sel = 1:c2
                imtemp = seq2masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.30*nrow(2)*ncol(2)) %select slices whose tumor part is more than 30%
                    arr{2} = [arr{2} sel];
                end
            end
        end
        if (isempty(arr{2}))
            for sel = 1:c2
                imtemp = seq2masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.20*nrow(2)*ncol(2)) %select slices whose tumor part is more than 30%
                    arr{2} = [arr{2} sel];
                end
            end
        end
        if (isempty(arr{2}))
            for sel = 1:c2
                imtemp = seq2masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.10*nrow(2)*ncol(2)) %select slices whose tumor part is more than 30%
                    arr{2} = [arr{2} sel];
                end
            end
        end
        if (isempty(arr{2}))
            for sel = 1:c2
                imtemp = seq2masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.04*nrow(2)*ncol(2)) %select slices whose tumor part is more than 30%
                    arr{2} = [arr{2} sel];
                end
            end
        end 
        if (isempty(arr{2}))
            for sel = 1:c2
                imtemp = seq2masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0*nrow(2)*ncol(2)) %select slices whose tumor part is more than 30%
                    arr{2} = [arr{2} sel];
                end
            end
        end 
        for sel = 1:c4
            imtemp = seq4masks{p}(:,:,sel);
            if(sum(imtemp(:)) > 0.50*nrow(3)*ncol(3)) %select slices whose tumor part is more than 30%
                arr{3} = [arr{3} sel];
            end
        end
        if (isempty(arr{3}))
            for sel = 1:c4
                imtemp = seq4masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.30*nrow(3)*ncol(3)) %select slices whose tumor part is more than 30%
                    arr{3} = [arr{3} sel];
                end
            end
        end
        if (isempty(arr{3}))
            for sel = 1:c4
                imtemp = seq4masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.20*nrow(3)*ncol(3)) %select slices whose tumor part is more than 30%
                    arr{3} = [arr{3} sel];
                end
            end
        end
        if (isempty(arr{3}))
            for sel = 1:c4
                imtemp = seq4masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.10*nrow(3)*ncol(3)) %select slices whose tumor part is more than 30%
                    arr{3} = [arr{3} sel];
                end
            end
        end
        if (isempty(arr{3}))
            for sel = 1:c4
                imtemp = seq4masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0.04*nrow(3)*ncol(3)) %select slices whose tumor part is more than 30%
                    arr{3} = [arr{3} sel];
                end
            end
        end 
        if (isempty(arr{3}))
            for sel = 1:c4
                imtemp = seq4masks{p}(:,:,sel);
                if(sum(imtemp(:)) > 0*nrow(3)*ncol(3)) %select slices whose tumor part is more than 30%
                    arr{3} = [arr{3} sel];
                end
            end
        end 
        numpat = [];
        numpat(1) = ceil(totalpat/length(arr{1})); %number of patches to sample per slice
        numpat(2) = ceil(totalpat/length(arr{2}));
        numpat(3) = ceil(totalpat/length(arr{3}));
        indx = 0;
        for i = 1:3
            arra = arr{i};
            indx = indx + 1;
            im = images{indx};
            npat = numpat(indx);
            tcount = 0; %total patch count
            tempallfeats = []; %temporary variable
            for slice = 1:length(arra)
                pcount = 0;
                while(1)
                    if(pcount > npat)
                        break;
                    end
                    row = randi(nrow(indx)-psize, 1);
                    col = randi(ncol(indx)-psize, 1);                    
                    patch = im(row:row+psize-1, col:col+psize-1, arra(slice));
                    pcount = pcount + 1; % patch count for present slice
                    tcount = tcount + 1;
                    tempallfeats(tcount,:) = patch(:);
                end
            end 
            allfeats{indx}{p} = tempallfeats(1:totalpat,:);
        end
    end
 
    %% Patch feature extraction
    n = length(finlabels);
    progressbar('CV repetitions', 'CV fold number');
    acc = []; auc = [];
    parfor repeat = 1:5 %repeat the cross-validation process multiple times
        cpart = cvpartition(finlabels,'KFold', 10);
        nfold = length(cpart.TrainSize); 
        for trial = 1:nfold
            trinds = find(cpart.training(trial)==1);
            tsinds = find(cpart.training(trial)==0);
            trfeats1 = []; tsfeats1 = [];
            trfeats2 = []; tsfeats2 = [];
            trfeats3 = []; tsfeats3 = [];
            for i = 1:length(trinds)
                trfeats1 = vertcat(trfeats1, allfeats{1}{trinds(i)});
                trfeats2 = vertcat(trfeats2, allfeats{2}{trinds(i)});
                trfeats3 = vertcat(trfeats3, allfeats{3}{trinds(i)});
            end
            for i = 1:length(tsinds)
                tsfeats1 = vertcat(tsfeats1, allfeats{1}{tsinds(i)});
                tsfeats2 = vertcat(tsfeats2, allfeats{2}{tsinds(i)});
                tsfeats3 = vertcat(tsfeats3, allfeats{3}{tsinds(i)});
            end
 
            %Dictionary Learning
            class11 = []; class01 = [];
            class12 = []; class02 = [];
            class13 = []; class03 = [];
            for i = 1:length(trinds)
                if(finlabels(trinds(i)) == 1)
                    class11 = [class11 ; trfeats1((i-1)*totalpat+1:i*totalpat,:)];
                    class12 = [class12 ; trfeats2((i-1)*totalpat+1:i*totalpat,:)];
                    class13 = [class13 ; trfeats3((i-1)*totalpat+1:i*totalpat,:)];
                else
                    class01 = [class01 ; trfeats1((i-1)*totalpat+1:i*totalpat,:)];
                    class02 = [class02 ; trfeats2((i-1)*totalpat+1:i*totalpat,:)];
                    class03 = [class03 ; trfeats3((i-1)*totalpat+1:i*totalpat,:)];
                end
            end
            Y1 = [class01 ; class11]';
            Y2 = [class02 ; class12]';
            Y3 = [class03 ; class13]';
            Yrange1 = [1, length(class01), length(class01)+length(class11)];
            Yts1 = tsfeats1';
            Yrange2 = [1, length(class02), length(class02)+length(class12)];
            Yts2 = tsfeats2';
            Yrange3 = [1, length(class03), length(class03)+length(class13)];
            Yts3 = tsfeats3';
            [~, pred1, ~,pred_tr1, ~] = runDLSI(Y1, Yrange1, Yts1, totalpat,1,0.1);
            [~, pred2, ~,pred_tr2, ~] = runDLSI(Y2, Yrange2, Yts2, totalpat,1,0.1);
            [~, pred3, ~,pred_tr3, ~] = runDLSI(Y3, Yrange3, Yts3, totalpat,1,0.1);
            
            ensemble_train = [pred_tr1 pred_tr2 pred_tr3];
 
            %fit ensemble
            model = fitcensemble(ensemble_train,finlabels(trinds),'method','AdaBoostM1','CategoricalPredictors','all','LearnRate',0.1);
            ensemble_test = [pred1 pred2 pred3]; 
            %make ensemble predictions
            [labels,score] = predict(model,ensemble_test)
            [~,~,~,AUC] = perfcurve(finlabels(tsinds),score(:,2),1);
            [pred1 pred2 pred3 labels finlabels(tsinds)]
            AUC
            auc = [auc AUC];
            progressbar([], trial/nfold);
        end
        progressbar(repeat/5, 0);
    end
    mean_aucs{labs} = auc;
end
end
