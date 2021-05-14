function [data_impute, W] = G2S3(data,varargin)
% Graph of smooth gene network signal inputed genes.
% 
% data format: cell(row) * genes(columns) array
%   'a' (default = 1)
%       bigger a -> bigger weights in gene networks
%   'b' (default = 1)
%       bigger b -> more dense W   
%   'mode'(default = 'logb')
%       'logb': with -log(degree_node) barrier to the learning of graph  
%        minimize_W sum(sum(W .* Z)) - a * sum(log(sum(W))) + b * ||W||_F^2/2 + c * ||W-W_0||_F^2/2
%
%
%       'l2p': add L2 penalty on the node degree
%        minimize_W sum(sum(W .* Z)) + a/2 * ||W||_F^2/2 
%                 + a/2 * ||sum(W)||_2^2 
%                 s.t. sum(sum(W)) == n
%   'tol'(default = 1e-4)
%       cut-off for neglectable weights   
%   'normalize' (default = 1)
%       whether or not normalize the cell library in the output data
%   'scale' (default = 1)
%       whether or not scale the gene expression level back
% output format: cell(row) * genes(columns) array

% set up default parameters   
a = 1; b = 1; type = 'logb';tol = 10^-4;self = 1; normalize = 1;scale = 1; t = 1;
% get input parameters

for i=1:length(varargin)
    % mode
    if(strcmp(varargin{i},'mode'))
       type = lower(varargin{i+1});
    end
    % a 
    if(strcmp(varargin{i},'a'))
       a = lower(varargin{i+1});
    end
    % b
    if(strcmp(varargin{i},'b'))
       b = lower(varargin{i+1});
    end
    % tol
    if(strcmp(varargin{i},'epsi'))
       tol = lower(varargin{i+1});
    end
    if(strcmp(varargin{i}, 'normalize'))
      normalize = lower(varargin{i+1});
    end

    if(strcmp(varargin{i}, 'self'))
      self = lower(varargin{i+1});
    end

    if(strcmp(varargin{i}, 'scale'))
      normalize = lower(varargin{i+1});
    end

    if(strcmp(varargin{i}, 't'))
      t = lower(varargin{i+1});
    end
  
end

[~, nc] = size(data);

% log transformation

my_normc = @(m)bsxfun(@rdivide,m,sqrt(sum(m.^2)));
M = my_normc(data);


function W = affinity2graph(Z)
    switch type
        case 'logb'
            [W, stat] = gsp_learn_graph_log_degrees(Z, a, b);
        case 'l2p'
            [W, stat] = gsp_learn_graph_l2_degrees(Z, a);

    W(W<tol) = 0;
    end
end


  Z = gsp_distanz(M).^2;
  W = affinity2graph(Z);

    
switch self
  case 1
    W = bsxfun(@rdivide, W, sum(W,1));
    W = W + eye(nc);
    disp('adding self prob')
  case 0
    W = W;
    disp('only use neighbor')
  end

W(W<0) = 0;

W = bsxfun(@rdivide, W, sum(W,1));
W = W^t;
data_impute = data * W;



if normalize ==1
  data_impute = bsxfun(@rdivide, data_impute, sum(data_impute,2))*mean(sum(data_impute,2));
  disp('normalized by library size')
else
  data_impute = data_impute;
  disp('unormalized')
end

if scale ==1
  data_impute = data_impute./sum(data_impute);
  data_impute = data_impute.*sum(data);
end

end
