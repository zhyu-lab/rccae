function segments = CloneCNA_segment_results(state_seq)
% 05/24/2022 by Zhenhua
% This function is used to generate segments from a sequence of assigned states

pre_state = [];
s_indx = [];
segments = [];
for i = 1:length(state_seq)
    if isempty(pre_state)
        pre_state = state_seq(i);
        s_indx = i;
    elseif state_seq(i) ~= pre_state
        segments = [segments;s_indx i-1 pre_state];
        pre_state = state_seq(i);
        s_indx = i;
    end
end

if s_indx <= length(state_seq)
    segments = [segments;s_indx length(state_seq) pre_state];
end

end