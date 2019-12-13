function outputscore = felsensteintree(nibbatree,normfreq,numleaves,mat_labels,leafnames,N,pointers)
possible = 1:numleaves;
phytree = nibbatree;
allscores = zeros(10000,phytree.NumNodes);
leafnodepunish = 0.1;


allpossparents = [];
possnodeidents = cell(phytree.NumNodes,1);
    clear i

while length(possible) < phytree.NumNodes
    for i = 1:length(pointers)
        if any(possible == pointers(i,1)) && any(possible == pointers(i,2)) %if the node is currently possible to be scored
            if pointers(i,1) <= numleaves && pointers(i,2) <= numleaves %if both branches are leaves
                node1 = pointers(i,1);
                node2 = pointers(i,2);
                barcode1 = mat_labels((find(leafnames==str2num(N{node1})==1)),:);
                barcode2 = mat_labels((find(leafnames==str2num(N{node2})==1)),:);
                similaritycoeff = (sum(barcode1==barcode2))^2; %Match the two leaf nodes
                possparents = intersect(allpossparent(barcode1),allpossparent(barcode2),'rows');
                allpossparents = unique(vertcat(allpossparents,possparents),'rows','stable');
                possnodeidents{i+numleaves} = possparents;
                for parent = 1:size(possparents,1)
                    ancestor = possparents(parent,:);
                    score = 1;
                    for j = 1:10
                        t = [1 0 0; normfreq(1,j) normfreq(2,j) normfreq(3,j); 0 0 1];
                        score = score * t(ancestor(j)+1,barcode1(j)+1)*t(ancestor(j)+1,barcode2(j)+1);
                    end

                    idx = find(ismember(allpossparents,ancestor,'rows')); 
                    allscores(idx,i+numleaves) = score*similaritycoeff;
                end
            elseif pointers(i,1) > numleaves && pointers(i,2) > numleaves %if both branches are internal nodes
                node1 = pointers(i,1);
                node2 = pointers(i,2);
                possnode1barcodes = possnodeidents{node1};
                possnode2barcodes = possnodeidents{node2};
                possnode1parents = [];
                possnode2parents = [];

                for x = 1:size(possnode1barcodes,1)
                    possnode1parents = vertcat(possnode1parents, allpossparent(possnode1barcodes(x,:)));
                end
                for x = 1:size(possnode2barcodes,1)
                    possnode2parents = vertcat(possnode2parents, allpossparent(possnode2barcodes(x,:)));
                end
                            
                possparents = intersect(possnode1parents,possnode2parents,'rows');
                allpossparents = unique(vertcat(allpossparents,possparents),'rows','stable');
                possnodeidents{i+numleaves} = possparents;
                for parent = 1:size(possparents,1)
                    parentbarcode = possparents(parent,:);
                    runsum1 = 0;
                    runsum2 = 0;
                    for child = 1:size(possnode1barcodes,1)
                        childbarcode = possnode1barcodes(child,:);
                        transcore = 1;
                        for j = 1:10
                            t = [1 0 0; normfreq(1,j) normfreq(2,j) normfreq(3,j); 0 0 1];
                            transcore = transcore * t(parentbarcode(j)+1,childbarcode(j)+1);
                        end
                        idx = find(ismember(allpossparents,childbarcode,'rows'));
                        runsum1 = runsum1 + allscores(idx,pointers(i,1))*transcore;
                    end
                    for child = 1:size(possnode2barcodes,1)
                        childbarcode = possnode2barcodes(child,:);
                        transcore = 1;
                        for j = 1:10
                            t = [1 0 0; normfreq(1,j) normfreq(2,j) normfreq(3,j); 0 0 1];
                            transcore = transcore * t(parentbarcode(j)+1,childbarcode(j)+1);
                        end
                        idx = find(ismember(allpossparents,childbarcode,'rows'));
                        runsum2 = runsum2 + allscores(idx,pointers(i,2))*transcore;
                    end
                    idx = find(ismember(allpossparents,parentbarcode,'rows'));
                    allscores(idx,i+numleaves) = runsum1*runsum2;
                end
            else %if one is leaf one is internal node
                nodenums = [pointers(i,1) pointers(i,2)];
                barcode1 = mat_labels(find(leafnames==str2num(N{min(nodenums)})),:); %the smaller index is the leaf
                possintnodeidents = possnodeidents{max(nodenums)};
                possintnodeparents = [];
                for x = 1:size(possintnodeidents,1);
                    possintnodeparents = vertcat(possintnodeparents, allpossparent(possintnodeidents(x,:)));
                end
                possparents = intersect(possintnodeparents,allpossparent(barcode1),'rows');
                allpossparents = unique(vertcat(allpossparents,possparents),'rows','stable');
                possnodeidents{i+numleaves} = possparents;
                possbignodebarcodes = possnodeidents{max(nodenums)};
                for parent = 1:size(possparents,1);
                    parentbarcode = possparents(parent,:);
                    leafscore = 1;
                    runsum = 0;
                    for j = 1:10
                        t = [1 0 0; normfreq(1,j) normfreq(2,j) normfreq(3,j); 0 0 1];
                        leafscore = leafscore * t(parentbarcode(j)+1,barcode1(j)+1);
                    end
                    for child = 1:size(possbignodebarcodes,1)
                        childbarcode = possbignodebarcodes(child,:);
                        transcore = 1;
                        for j = 1:10
                            t = [1 0 0; normfreq(1,j) normfreq(2,j) normfreq(3,j); 0 0 1];
                            transcore = transcore * t(parentbarcode(j)+1,childbarcode(j)+1);
                        end
                        idx = find(ismember(allpossparents,childbarcode,'rows'));
                        runsum = runsum + allscores(idx,max(nodenums))*transcore;
                    end
%                     if numleavenode > 3
%                         leafnodepunish = 0;
%                     end
                    idx = find(ismember(allpossparents,parentbarcode,'rows'));
                    allscores(idx,i+numleaves)= leafscore*runsum*leafnodepunish;
                end
            end
            possible = [possible i+numleaves];
            
        end
    end
end

rownum = find(ismember(allpossparents,[1 1 1 1 1 1 1 1 1 1],'rows'));
colsum = sum(allscores,1);
outputscore = colsum(end);
%outputscore = allscores(rownum,end);
end