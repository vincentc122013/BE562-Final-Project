function possparents = allpossparent(barcode)
possparents = zeros(1,length(barcode));

for j = 1:length(barcode)
    if barcode(j) == 1
        possparents(:,j) = 1;
    else
        possparents = vertcat(possparents, possparents);
        possparents(1:size(possparents,1)/2,j) = barcode(j);
        possparents(size(possparents,1)/2 + 1:size(possparents,1),j) = 1;
    end
end