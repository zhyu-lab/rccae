function results = SCHMM_eval_pdf_lrc(data,mu,sigma)

if size(data,1) > size(data,2) %Nx1->1xN
    data = data';
end

results = normpdf(data,mu,sigma);

end