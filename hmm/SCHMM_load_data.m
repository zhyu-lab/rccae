function [data_lrc_all,data_chr_all,data_bin_all,bin_size] = SCHMM_load_data(dataFile)
% function [data_lrc_clone,data_chr_all,data_bin_all,barcords,bin_size,labels_n] = CloneCNA_load_data(matFile,labelFile)

% load(matFile,'data_lrc_all','data_chr_all','data_bin_all','barcords','bin_size');
fid = fopen(dataFile, 'r');
line = fgetl(fid);
bin_size = str2double(line);
line = fgetl(fid);
fields = regexp(line,',','split');
data_chr_all = str2double(fields);
line = fgetl(fid);
fields = regexp(line,',','split');
data_bin_all = str2double(fields);
results = textscan(fid,repmat('%f',1,length(data_chr_all)),'Delimiter',',');
data_lrc_all = cell2mat(results);
clear results;
fclose(fid);

data_rc_all = 2.^data_lrc_all;
for i = 1:size(data_rc_all,1)
    data_rc_all(i,:) = data_rc_all(i,:)/(median(data_rc_all(i,:))+eps);
end
data_lrc_all = log2(data_rc_all+eps);

chromosomes = reshape(unique(data_chr_all),1,[]);
sorted_indxs = [];
for i = 1:length(chromosomes)
    indxs = find(data_chr_all == chromosomes(i));
    sorted_indxs = [sorted_indxs indxs];
end

data_lrc_all = data_lrc_all(:,sorted_indxs);
data_chr_all = data_chr_all(sorted_indxs);
data_bin_all = data_bin_all(sorted_indxs);

end