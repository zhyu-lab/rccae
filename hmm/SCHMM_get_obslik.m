function [obslik,condi_probs_fluct] = SCHMM_get_obslik(data_lrc,o,sigma,depend_table)

N = length(data_lrc); %number of data points

tv_S = depend_table(:,2) == 1;
Y = depend_table(tv_S,3); %vector of copy numbers of different entries
mu_l = log2(Y/2)+o;

S = sum(tv_S);
obslik = zeros(S,N);
condi_probs_fluct = zeros(S,N);

fluct_prob = 1e-5;

for i = 1:length(Y)
    obslik_lrc = SCHMM_eval_pdf_lrc(data_lrc,mu_l(i),sigma);
    if Y(i) == 2
        obslik_lrc = 1.06*obslik_lrc;
    end
    if Y(i) == 1
        obslik_lrc = 1.2*obslik_lrc;
    end
    obslik(i,:) = (1-fluct_prob)*obslik_lrc+fluct_prob/6;
    condi_probs_fluct(i,:) = (fluct_prob/6)./obslik(i,:);
end

end

