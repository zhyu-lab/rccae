function [p_states,aCN,segments] = SCHMM_process_results(depend_table)
% 05/24/2022 by Zhenhua
%-----------------------------------------------------
%------overall information of the clone------
%p_states: proportions of all hidden states
%num_loci: total number of loci investigated
%aCN: averaged copy number
%segments: copy number segmentation results

global gamma_sep

tv_S = depend_table(:,2)==1;
Y = depend_table(tv_S,3)'; % copy number of different entries

%initialize output parameters
segments = [];

%initialize intermediate variables
exp_num_states = [];
pos_dist = [];

for i = 1:length(gamma_sep) %for the ith chromosome
    post_probs = gamma_sep{i};
    
    %---handle p_states and num_loci---
    if isempty(exp_num_states) %initialization
        exp_num_states = zeros(size(post_probs,1),1);
    end
    exp_num_states = exp_num_states+sum(post_probs,2);

    %---handle MAP states---
    %output predicted MAP states
    [temp,MAP_state] = max(post_probs,[],1);

    results = SCHMM_segment_results(MAP_state);
    segments = [segments; ones(size(results,1),1)*i results];    
    pos_dist = [pos_dist; (results(:,2)-results(:,1)+1)]; 

    clear results;

end
pos_dist = pos_dist';

%---handle p_states---
p_states = zeros(length(Y),1);
for i = 1:length(Y)
    tv = segments(:,4) == i;
    if sum(tv) > 0
        p_states(i) = sum(pos_dist(tv))/sum(pos_dist);
    end
end

%---handle aCN---
aCN = Y*p_states;

end