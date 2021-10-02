function [ labels, mdata ] = quantizeNL_IRK( y, nclust )
% y is input data
% nclust is the number of clusters or quantization levels. If you want to
% quantize data with b bits then nclust <= 2^b
% labels are the quantization intervals in [1, nclust] range assigned to y
% mdata is thhe mean value of each cluster

y = double(y);
y0 = y;
y = reshape(y, 1, numel(y));
[y, idxY] = sort(y);

edges = [0,numel(y)];
errors = sum((y-mean(y)).^2);

s_data = cumsum(y);
ss_data = cumsum(y.^2);

for i = 1 : nclust-1
    [~,idx] = max(errors);
    k = edges(idx); n = edges(idx+1)-k;
    sn = s_data(k+n); if(k>=1) sn = sn - s_data(k); end
    ssn = ss_data(k+n); if(k>=1) ssn = ssn - ss_data(k); end
    d = 2; m = floor(n/d);
    while(1)
        sm = s_data(k+m); if(k>=1) sm = sm - s_data(k); end
        ssm = ss_data(k+m); if(k>=1) ssm = ssm - ss_data(k); end
        e1 = ssm-sm^2/m;
        e2 = ssn - ssm - (sn - sm)^2/(n-m);
        d = 2 * d;
        if(abs(e1-e2) < 0.001 || d >= n)
            edges = [edges(1:idx),k+m,edges(idx+1:end)];
            errors = [errors(1:idx-1),e1,e2,errors(idx+1:end)];
            break;
        else
            if(e1 > e2) m = m-floor(n/d); elseif(e1 < e2) m = m+floor(n/d); end
        end
    end
end

labels = [];
j = 1;
for i=1:numel(edges)-1
    labels = [labels,j*ones(1,edges(i+1)-edges(i))];
    j = j+1;
end
labels(idxY) = labels;
labels = reshape(labels, size(y0));

mdata = zeros(1, nclust);
for i = 1:nclust
    j = (labels == i);
    mdata(i) = mean(y0(j));
end
end

