load wspace.mat;
%% Select patients having a particular MR sequence available
inds = {[],[],[],[]}; %Flair, T1, T1C, T2
for j = 1:4
    for i = 1:length(allimagesROI{1})
        if(~isempty(allimagesROI{j}{i}))
            inds{j} = [inds{j} , i];
        end
    end
end

seq = 1; %Select sequence: 1.Flair, 2.T1, 3.T1C, 4.T2
imcells = {allimagesROI{seq}{inds{seq}}};
masks = {allimagesmasks{seq}{inds{seq}}};

for i = 1:length(imcells)
    minimum = prctile(imcells{i}(:),1);    %( 3)
    maximum = prctile(imcells{i}(:),99);
    imcells{i} = (imcells{i} - minimum)/(maximum - minimum);
end

labels = labels(inds{seq});
glabels = glabels(inds{seq});
clabels = clabels(inds{seq});
hlabels = hlabels(inds{seq});
slabels = slabels(inds{seq});
allpws = allpws(inds{seq},:);

%{
seq = 2;
masks = {allimagesmasks{seq}{inds{seq}}};
arr = [];
for i = 1:length(masks)
    arr = [arr min(size(masks{i},1) , size(masks{i},2))];
end
min(arr)  
%}



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
%% All patients training
labelsall = [labels, glabels, clabels];
for labs = 2
    finlabels = labelsall(:,labs);
    totalpat = 500; %number of patches per patient
    allfeats = {}; %contains all the image patches 
    psize = 20; %patch size: [psize x psize]
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
        [arrsort, inds] = sort(arr, 'descend');
        npat = ceil(totalpat/nslices); %number of patches to sample per slice
        inds = inds(1:nslices)
        tcount = 0; %total patch count
        tempallfeats = []; %temporary variable
        for slice = inds
            pcount = 0;
            while(1)
                if(pcount> npat)
                    break;
                end
                row = randi(nrow-psize, 1);
                col = randi(ncol-psize, 1);
                patch = im(row:row+psize-1, col:col+psize-1, slice);
                pcount = pcount + 1; % patch count for present slice
                tcount = tcount + 1;
                tempallfeats(tcount,:) = patch(:);
            end
        end 
        allfeats{p} = tempallfeats(1:totalpat,:);
    end

    %% Model training and cross validation
    n = length(finlabels);
    progressbar('CV repetitions', 'CV fold number');
    for repeat = 1:5 %repeat the cross-validation process multiple times
        acc = []; auc = [];
        cpart = cvpartition(finlabels,'KFold', 5);
        nfold = length(cpart.TrainSize); 
        for trial = 1:nfold
            trinds = find(cpart.training(trial)==1);
            tsinds = find(cpart.training(trial)==0);
            trfeats = []; tsfeats = [];
            for i = 1:length(trinds)
                trfeats = vertcat(trfeats, allfeats{trinds(i)});
            end
            for i = 1:length(tsinds)
                tsfeats = vertcat(tsfeats, allfeats{tsinds(i)});
            end
            
            %Refer to bag-of-words algorithm in computer vision
            %K-Means
            opt = statset('UseParallel',1);
            N = 10; %User-defined parameter [choose best value !!]
            tic;
            [idx,C,sumd,D] = kmeans(trfeats, N, 'Maxiter', 1000,'options', opt);
            toc;

            %Binning for training set
            trbinvec = []; %For each patient, calculate a histogram (feature vector)
            for i = 1:length(trinds)
                lookvec = (i-1)*totalpat + 1 : i*totalpat;
                [a1, a2] = hist(idx( lookvec ), 1:N);
                trbinvec(i,:) = a1;
            end

            %{
            % Copy-paste this chunk to View cluster centres
            for i = 1:size(C)
                i
                patch = reshape(C(i,:), [psize, psize]);
                imtemp = (patch - min(patch(:)))/(max(patch(:)) - min(patch(:)));
                imshow(imtemp, 'InitialMagnification', 500);
                pause;
            end
            %}

            %Test set
            assgn = [];
            for i = 1:length(tsfeats)
                %i
                a = pdist2(tsfeats(i,:), C); %Distances of test points from cluster centers
                [a1, a2] = min(a);
                assgn = [assgn a2];
            end
            % Test bins
            tsbinvec = [];
            for i = 1:length(tsinds)
                lookvec = (i-1)*totalpat + 1 : i*totalpat;
                [a1, a2] = hist(assgn( lookvec ), 1:N);
                tsbinvec(i,:) = a1;
            end
            
            %SVM model training [can choose any classifier]
            Xtr = trbinvec;
            ytr = glabels(trinds);
            mdl = fitcsvm(Xtr, ytr, 'KernelFunction', 'poly');

            Xts = tsbinvec;
            yts = glabels(tsinds);
            
            [pred score cost] = predict(mdl, Xts);
            [pred finlabels(tsinds) score(:,2)]
            acc = [acc sum(pred == yts)/length(pred) ];
            [x1,x2,T,AUC] = perfcurve(yts,score(:,2),1);
            auc = [auc AUC];
            AUC
            progressbar([], trial/nfold);
        end
        progressbar(repeat/5, 0);
        mean(auc)
    end
end
