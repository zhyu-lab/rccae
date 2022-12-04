function SCHMM_plot(inputFile,resultDir)

plotDir = [resultDir '/plots'];
mkdir(plotDir);

% load data
[data_lrc_all,data_chr_all,data_bin_all,bin_size] = SCHMM_load_data(inputFile);
data_spos_all = (data_bin_all-1)*bin_size+1;
data_epos_all = data_bin_all*bin_size;

num_cell = size(data_lrc_all,1);

fid = fopen([resultDir '/paras.csv'],'r');
if fid == -1
    error(['Can not open file: ' resultDir '/paras.csv']);
end
results = textscan(fid,'%f%*f','headerlines',1,'delimiter',',');
fclose(fid);
o_all = results{1};

fid = fopen([resultDir '/ploidy.csv'],'r');
if fid == -1
    error(['Can not open file: ' resultDir '/ploidy.csv']);
end
line = fgetl(fid);
fields = regexp(line,',','split');
ploidy_all = str2double(fields);
fclose(fid);

fid = fopen([resultDir '/segments.csv'],'r');
if fid == -1
    error(['Can not open file: ' resultDir '/segments.csv']);
end
results = textscan(fid,'%f%f%f%f%f','headerlines',1,'delimiter',',');
fclose(fid);
cell_seg_all = results{1};
chr_seg_all = results{2};
spos_seg_all = results{3};
epos_seg_all = results{4};
cn_seg_all = results{5};

% for k = 833:num_cell
for k = 1:num_cell
    tv = cell_seg_all == k;
    chr_seg = chr_seg_all(tv);
    pstart_seg = spos_seg_all(tv);
    pend_seg = epos_seg_all(tv);
    cn_seg = cn_seg_all(tv);
    
    o = o_all(k);
    ploidy = ploidy_all(k);

    h = figure(1);
    set(h,'visible','off');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 4])
    FontSize = 11;

    %--------------- plot figures ---------------------%
    lcr_colors = [0.5 0.5 0.5;
                 0 0.9 0;
                 0 0 0.9;
                 0.9 0 0];

    chromosomes = intersect(unique(data_chr_all),1:22);
    max_pos = zeros(1,length(chromosomes));
    
    data_lrc_cell = data_lrc_all(k,:);
%     max_lrc = max(data_lrc_cell);
%     min_lcr = min(data_lrc_cell);

    for i = 1:length(chromosomes)
        tv1 = data_chr_all == chromosomes(i);
        max_pos(i) = max(data_epos_all(tv1));
    end

    ratio = max_pos/sum(max_pos);
    xtick = cumsum([0 ratio(1:end-1)])+ratio/2;

    clf;
    line_style = 'k-';
    LineWidth = 1.0;
    MarkerSize = 4;
    subplot(2,1,1);
    hold on
    set(gca,'YGrid','on');
    set(gca,'FontSize',FontSize);
    set(gca,'YTick',[0:1:7],'Box','on');
    set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
    ylabel('Copy number');
    pre_x = 0;
    chr_epos = zeros(length(chromosomes),1);
    for i = 1:length(chromosomes)
        tv = data_chr_all == chromosomes(i);
        data_spos = data_spos_all(tv);
        data_epos = data_epos_all(tv);
        x = data_epos*ratio(i)/max_pos(i)+pre_x;
        indx1 = find(chr_seg == chromosomes(i));
        for j = reshape(indx1,1,[])
            CN = cn_seg(j);
            indx = find(data_spos >= pstart_seg(j) & data_epos <= pend_seg(j));
            if isempty(indx)
                continue;
            end
            if CN > 7
                CN = 7;
            end
            plot([x(indx(1)) x(indx(end))],[CN CN], 'r-', 'LineWidth',2.0);
        end
        chr_epos(i) = max(x);
        pre_x = pre_x+ratio(i);  
    end
    for i = 1:length(chromosomes)-1
        plot([chr_epos(i) chr_epos(i)], [-0.1 7.1], line_style, 'LineWidth',LineWidth)
    end
    set(gca,'XTick',[])
    axis([0 1 -0.1 7.1]);
    tmp = ['Cell ' num2str(k) ', ploidy = ' num2str(ploidy)];
    title(tmp);

    subplot(2,1,2);
    set(gca,'FontSize',FontSize);
    hold on
    pre_x = 0;
    for i = 1:length(chromosomes)
        tv = data_chr_all == chromosomes(i);
        data_lcr = data_lrc_cell(tv);
        data_spos = data_spos_all(tv);
        data_epos = data_epos_all(tv);
        x = data_epos*ratio(i)/max_pos(i)+pre_x;
        indx1 = find(chr_seg == chromosomes(i));
        for j = reshape(indx1,1,[])
            CN = cn_seg(j);
            tv = data_spos >= pstart_seg(j) & data_epos <= pend_seg(j);
            if sum(tv) == 0
                continue;
            end
            if CN < 1
                c = 1;
            else
                c = CN+1;
            end
            if c > 4
                c = 4;
            end
            plot(x(tv),data_lcr(tv),'.','MarkerSize',MarkerSize, 'Color', lcr_colors(c,:));
        end
        % plot expected LCR mean values
        for j = reshape(indx1,1,[])
            CN = cn_seg(j);
            if CN == 0
                CN = 0.001;
            end
            lrc_mean = log2(CN/2)+o;
            indx = find(data_spos >= pstart_seg(j) & data_epos <= pend_seg(j));
            if isempty(indx)
                continue;
            end
            plot([x(indx(1)) x(indx(end))],[lrc_mean lrc_mean],'k-','LineWidth',1.5);
        end
        pre_x = pre_x+ratio(i);  
    end
    for i = 1:length(chromosomes)-1
        plot([chr_epos(i) chr_epos(i)], [-3 3], line_style, 'LineWidth',LineWidth)
    end
    ylabel('LRC');
    set(gca,'Box','on')
    axis([0 1 -2 2])
    xlabel('Chromosome');
    set(gca,'XTick',xtick);
    set(gca,'XTickLabel',mat2cell(chromosomes,1,length(chromosomes)));

    %save figure
    figpath = [plotDir '/cell_' num2str(k) '.png'];
    eval(['print -dpng -r100 ' figpath ])
end

end