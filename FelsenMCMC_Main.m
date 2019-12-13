clear; clc; close all;
%Read all subtrain file
allbarcodes = readmatrix('allBarcodes.txt');

for ifile = 1:50
   filename = ['sub1_train_',num2str(ifile)];
   trainfile = [filename '.txt'];
   groundfile = [filename '.nwk'];
   outputfile = [num2str(ifile) 'outtree.nwk'];
%% initialize
allpointer = [];
allscore = [];
train1 = importdata(trainfile); % Data for train 1
train1_GT = importdata(groundfile); % Ground truth for train 1
allN = {};
train1_GTtree = phytreeread('sub1_train_1.nwk'); % Phylogenetic tree from ground truth of train 1
%% mat_labels is a matrix that holds cell state values
initialcellname = train1.data(:,1);
train1.data(:,1) = 1:length(train1.data);
leafnames = train1.data(:,1);
mat_labels = zeros(length(train1.data(:,2)),10);
valindex = [];
for i = 1:length(train1.data(:,2))
    val = num2str(train1.data(i,2));
    indexnum = length(val);
    extrazeros = num2str(zeros(1,10-indexnum));
    val = [extrazeros val];
    val = val((~isspace(val)));
    valindex = [valindex;val];
    for j = 1:10 
        mat_labels(i,j) = val(j);
    end
end
mat_labels = mat_labels-48;
for ii = 1:length(train1.data(:,1))
    cella{ii,1} = num2str(train1.data(ii,1));
    cella{ii,2} = valindex(ii,:);
end
fields = {'Header', 'Sequence'};
cellstruct = cell2struct(cella, fields, 2);

numleaves = size(mat_labels,1);
%% Implement UPGMA
dist = seqpdist(cellstruct,'method','jukes-cantor','indels','pair');
orig_tree =  seqlinkage(dist,'average',cellstruct);
%%
freq=zeros(3,10);
treei = get(orig_tree);
for i = 1:10
    freq(2,i) = sum((allbarcodes(:,i) == 1)); %ones
    freq(3,i) = sum((allbarcodes(:,i) == 2)); %twos
    freq(1,i) = sum((allbarcodes(:,i) == 0)); %zeros
end
normfreq = freq/1029;

%%For loop starts here to iterate and file optimal tree
for runnum = 1:200
    
    %% Metropolis MCMC
    pointers = treei.Pointers;
    newp = pointers;
    N = treei.NodeNames;
    for numN = 1:length(cellstruct)
        N{numN}=char(regexp(N{numN}, '\d+', 'match'));
    end
    oldtree = phytree(pointers);
    
    
    
    nodechildren = newp(randi([1,length(newp)-1]),:); %picks the row containing the children of one of the internal nodes, but not the root node
    childrenrowidx = find(pointers==nodechildren,1); %find the row containing the children of your internal node
    nodenum = numleaves + childrenrowidx; %find the number of your internal node
    noderowidx = find(any(pointers==nodenum,2)); %find the row containing your node and its sister
    
    sistercolidx = find(pointers(noderowidx,:) ~= nodenum);
 
    picked = randi([1 2]);
    newp(childrenrowidx,picked)=pointers(noderowidx,sistercolidx);
    newp(noderowidx,sistercolidx)=pointers(childrenrowidx,picked);
    
    
    %Code to fix tree
    leavesize = size(mat_labels,1);
for iii=1:length(newp)
    threshold = leavesize+iii;
    if max(newp(iii,:)) >= threshold
        maxval = max(newp(iii,:));
        bruhval = newp(iii,:);
        newp(iii,:) = newp(iii+1,:);
        newp(iii+1,:) = bruhval;
        [r1,c1]=find(newp==maxval);
        [r2,c2]=find(newp==threshold);
        newp(r1,c1) = threshold;
        newp(r2,c2) = maxval;
        disp('tree pointer swapped to fix order');
    end
end
    newtree = phytree(newp,N);
    getnewtree = get(newtree);
    Nn = getnewtree.NodeNames;
        for numN = 1:length(cellstruct)
        Nn{numN}=char(regexp(Nn{numN}, '\d+', 'match'));
    end
    
    %% Felsenstein's On the Initial Tree
    if runnum > 1
    oldscore = allscore(runnum-1);
    else
    oldscore = felsensteintree(treei,normfreq,numleaves,mat_labels,leafnames,N,pointers);
    end
    disp('Old score:')
    disp(oldscore)
    clear i
    %% Felsenstein's On the New Tree
    newscore = felsensteintree(getnewtree,normfreq,numleaves,mat_labels,leafnames,Nn,newp);
    clear i
    disp('New score:')
    disp(newscore)
    
    
    logscore = log(newscore/oldscore);
    
  
    if logscore >= 0
        allpointer(:,:,runnum) = newp;
        treei = getnewtree;
        allscore(runnum) = newscore;
        disp('higher score obtained on new tree, swapped')
        
    else
        changeppbly = 1/(1+exp(-logscore));
        changebenchmark = rand;
        if changeppbly >= changebenchmark*0.5
            treei = getnewtree;
            allpointer(:,:,runnum) = newp;
            allscore(runnum) = newscore;
            disp('no higher score obtained on new tree, swapped due to probabilty')
        else
            allpointer(:,:,runnum) = pointers;
            allscore(runnum) = oldscore;
            disp('No swap')
        end
    end
disp('Iteration number:')
disp(runnum)

end

endpointer = allpointer(:,:,end);
endN = N;
for i = 1:length(cellstruct)
    endN{i} = num2str(initialcellname(i));
end
%a = getnewickstr(phytree(nibba))
for ii = 1:length(cellstruct)
cella{ii,1} = (endN{ii});
end

newN = endN;
for i = 1:numleaves
    cellnum = newN{i};
    idx = find(strcmp(cella(:,1), cellnum));
    newN{i} = strcat(cellnum,'_',cellstruct(idx).Sequence);
end
n = getnewickstr(phytree(newp,newN));
%pp = best p
%n = getnewickstr(phytree(pp,newN))

dlmwrite(outputfile,n,'delimiter','');

clearvars -except allbarcodes ifile
clc
end