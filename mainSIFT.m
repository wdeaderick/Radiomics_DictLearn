function [mean_aucs] = mainSIFT(method, seq)
%The method argument should be passed 1 for DLSI, 2 for COPAR, 3 for FDDL
%The seq argument should be passed 1 for the Flair modality, 2 for the T1 modality, 3 for the 
%T1C modality, 4 for the T2 modality
%Returns a 1x3 vectors "mean_aucs" containing the mean AUCs for IDH Status, Grade, and Codeletion, respectively.

load wspace.mat
mean_aucs = [];%Preallocate return vector
 
%{ Select patients having a particular MR sequence available
inds = {[],[],[],[]}; %Flair, T1, T1C, T2
for j = 1:4
    for i = 1:length(allimagesROI{1})
        if(~isempty(allimagesROI{j}{i}))
            inds{j} = [inds{j} , i];
        end
    end
end
 
%{
seq = 1; %Select sequence: 1.Flair, 2.T1, 3.T1C, 4.T2
masks = {allimagesmasks{seq}{inds{seq}}};
 
for i = 1:length(imcells)
    minimum = prctile(imcells{i}(:),1);    %( 3)
    maximum = prctile(imcells{i}(:),99);
    imcells{i} = (imcells{i} - minimum)/(maximum - minimum);
end
 
%}

masks = {allimagesmasks{seq}{inds{seq}}};
imcells = imcells_seq4; %Select sequence: 1.Flair, 2.T1, 3.T1C, 4.T2
 
labels = labels(inds{seq});
glabels = glabels(inds{seq});
clabels = clabels(inds{seq});
hlabels = hlabels(inds{seq});
slabels = slabels(inds{seq});
allpws = allpws(inds{seq},:);  
 
%% For survival - ignore label:2
%{
ind = [];
for i = 1:length(slabels)
    if(slabels(i) ~= 2)
        ind = [ind i];
    end
end
imcells = imcells(ind);
slabels = slabels(ind);
[a,b] = grp2idx(slabels);
a = a-1;
slabels = a;
finlabels = slabels;
%}
 
%}
%% All patients training
labelsall = [labels, glabels, clabels];
for labs = 1:3
    finlabels = labelsall(:,labs);
    totalpat = 140; %number of patches per patient
    allfeats = {}; %contains all the image patches 
    psize = 16; %patch size: [psize x psize]
    for p = 1:length(finlabels)
        im = double(imcells{p});
        nrow = size(im,1);
        ncol = size(im,2);
        c = size(im, 3);
        arr = [];
        for sel = 1:c
            imtemp = masks{p}(:,:,sel);
            arr = [arr sum(imtemp(:))];
        end
        nslices = min(length(arr), 5); %Choose best 5 slices with max tumor area
        [~, inds] = sort(arr, 'descend');
        npat = ceil(totalpat/nslices); %number of patches to sample per slice
        inds = inds(1:nslices);
        tcount = 1; %total patch count
        tempallfeats = []; %temporary variable
        for slice = inds
            [~,d] = vl_dsift(single(im(:,:,slice)), 'size', 4); %SIFT descriptor
            pcount = 1;
            while(1)
                if(pcount > min(size(d,2),npat))
                    break;
                end
                tempallfeats(tcount,:) = d(:,pcount);
                pcount = pcount + 1; % patch count for present slice
                tcount = tcount + 1;
            end
        end
        allfeats{p} = tempallfeats(1:end,:);
    end
    
    %imshow(im(:,:,slice))    
    %h1 = vl_plotframe(f(:,:)) ;
    %h2 = vl_plotframe(f(:,:)) ;
    %set(h1,'color','k','linewidth',3) ;
    %set(h2,'color','y','linewidth',2) ;
    %% Patch feature extraction
    n = length(finlabels);
    acc = []; auc = [];
    progressbar('CV repetitions', 'CV fold number');
    parfor repeat = 1:5 %repeat the cross-validation process multiple times
        cpart = cvpartition(finlabels,'KFold', 10);
        nfold = length(cpart.TrainSize); 
        for trial = 1:nfold
            trinds = find(cpart.training(trial)==1);
            tsinds = find(cpart.training(trial)==0);
            trfeats = []; tsfeats = []; trlengths = []; tslengths = [];
            for i = 1:length(trinds)
                trfeats = vertcat(trfeats, allfeats{trinds(i)});
                trlengths(i) = size(allfeats{trinds(i)},1);
            end
            for i = 1:length(tsinds)
                tsfeats = vertcat(tsfeats, allfeats{tsinds(i)});
                tslengths(i) = size(allfeats{tsinds(i)},1);
            end
            %Dictionary Learning
            class1 = []; class0 = [];
            bookmark = 0;
            for i = 1:length(trinds)
                if(finlabels(trinds(i)) == 1)
                    class1 = [class1 ; trfeats(bookmark + 1:bookmark + trlengths(i),:)];
                else
                    class0 = [class0 ; trfeats(bookmark + 1:bookmark + trlengths(i),:)];
                end
                bookmark = bookmark + trlengths(i);
            end
            Y = [class0 ; class1]';
            Yrange = [1, length(class0), length(class0)+length(class1)];
            Yts = tsfeats';
            if(method == 1)
                [~, ~, f1s, ~] = siftrunDLSI(Y, Yrange, Yts, tslengths);
            elseif(method == 2)    
                [~, ~, f1s] = siftrunCOPAR(Y, Yrange, Yts, tslengths);
            elseif(method == 3)    
                [~, ~, f1s] = siftrunFDDL(Y, Yrange, Yts, tslengths);
            else
                disp('Invalid Argument')
                break
            end
            %acc = [acc sum(pred == finlabels(tsinds))/length(pred) ];
            %AUC for test set
            [~,~,~,AUC] = perfcurve(finlabels(tsinds),f1s./(tslengths'),1);
            AUC
            auc = [auc AUC];
            progressbar([], trial/nfold);
        end
        progressbar(repeat/5, 0);
    end
    mean_aucs(labs) = mean(auc);
end
end
