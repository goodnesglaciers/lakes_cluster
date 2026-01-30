function bins = discretize(x, edges, varargin)
% discretize  Group numeric data into bins or categories.
%    BINS = discretize(X,EDGES) returns the indices of the bins that the 
%    elements of X fall into.  EDGES is a numeric vector that contains bin 
%    edges in monotonically increasing order. An element X(i) falls into 
%    the j-th bin if EDGES(j) <= X(i) < EDGES(j+1), for 1 <= j < N where 
%    N is the number of bins and length(EDGES) = N+1. The last bin includes 
%    the right edge such that it contains EDGES(N) <= X(i) <= EDGES(N+1). For 
%    out-of-range values where X(i) < EDGES(1) or X(i) > EDGES(N+1) or 
%    isnan(X(i)), BINS(i) returns NaN.
%     
%    BINS = discretize(X,EDGES,VALUES) returns the corresponding element in 
%    VALUES, rather than the bin number. For example, if X(1) falls into 
%    bin 5 then BINS(1) would be VALUES(5) rather than 5. VALUES must be 
%    a vector with length equal to the number of bins. Out-of-range inputs 
%    return NaN. 
%     
%    C = discretize(X,EDGES,'categorical') creates a categorical array from 
%    the binned result of X. C will have the same number of categories 
%    as the number of bins. The category names will be in the form of 
%    "[A,B)", or "[A,B]" for the last bin, where A and B are consecutive 
%    values from EDGES. Out-of-range values will be undefined in C.
%     
%    C = discretize(X,EDGES,'categorical',CATEGORYNAMES) creates a categorical 
%    array from the binned result of X, and names the categories in C using 
%    CATEGORYNAMES.  CATEGORYNAMES is a cell array of strings and has
%    length equal to the number of bins.  
% 
%    BINS = discretize(...,'IncludedEdge',SIDE) and 
%    C = discretize(...,'IncludedEdge',SIDE) 
%    specify which bin edge is included in the bins. SIDE can be:
%        'left'    Each bin includes the left bin edge, except for the last bin
%                  which includes both bin edges. This is the default.
%        'right'   Each bin includes the right bin edge, except for the first 
%                  bin which includes both bin edges.
%    If SIDE is 'right', an element X(i) falls into the j-th bin if
%    EDGES(j) < X(i) <= EDGES(j+1), for 1 < j <= N where N is the number
%    of bins. The first bin includes the left edge such that it contains 
%    EDGES(1) <= X(i) <= EDGES(2).
%    
%    Class support for inputs X, EDGES:
%       float: double, single
%       integers: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%       logical
% 
%    Example:
%       x = rand(10,1);
%       bins = discretize(x,0:0.25:1);  % group into 4 bins 
%
%    See also HISTCOUNTS, HISTOGRAM, CATEGORICAL

%   Copyright 1984-2015 The MathWorks, Inc.

nin = nargin;

funcname = mfilename();
validateattributes(x, {'numeric','logical'}, {'real'}, funcname, 'x', 1)
validateattributes(edges, {'numeric','logical'},{'vector', 'real', ...
    'nondecreasing'}, funcname, 'edges', 2)
if length(edges) < 2
    error(message('MATLAB:discretize:EmptyOrScalarEdges'));
end

% make sure edges are non-sparse, handle subclass of builtin class
if isobject(x)
    x = matlab.internal.language.castToBuiltinSuperclass(x);
end
if isobject(edges)
    edges = matlab.internal.language.castToBuiltinSuperclass(edges);
end
edges = full(edges);
nbins = length(edges)-1;
    
persistent p p2;
    
if nin > 2 && isrow(varargin{1}) && ~iscell(varargin{1}) ...
        && strncmpi(varargin{1},'categorical',max(length(varargin{1}),1))
    % create categorical output
    if nin > 3
        if isempty(p)
            p = inputParser;
            addOptional(p, 'categorynames', NaN, @(x) iscellstr(x) && isvector(x))
            addParameter(p, 'IncludedEdge', 'left', ...
                @(x) validateattributes(x,{'char'},{}))
        end
        parse(p,varargin{2:end})
        catnames = p.Results.categorynames;
        catnames_provided = iscell(catnames);
        if catnames_provided && length(catnames) ~= nbins
            error(message('MATLAB:discretize:CategoryNamesInvalidSize',nbins));
        end
        
        right = strcmp(validatestring(...
            p.Results.IncludedEdge,{'left','right'}),'right');
    else
        catnames_provided = false;
        right = false;
    end
    
    if ~catnames_provided
        catnames = gencatnames(edges,right);
    end
    
    bins = discretizemex(x, edges, right);
    
    bins = categorical(bins, 1:nbins, catnames, 'Ordinal', true);
else
    % create numerical output
    if nin > 2
        if isempty(p2)
            p2 = inputParser;
            addOptional(p2, 'values', [], @(x) validateattributes(x, {'numeric'}, ...
                {'vector','nonempty'}))
            addParameter(p2, 'IncludedEdge', 'left', ...
                @(x) validateattributes(x,{'char'},{}))
        end
        parse(p2,varargin{:})
        values = p2.Results.values;
        values_provided = ~isempty(values);
        if values_provided && length(values) ~= nbins
            error(message('MATLAB:discretize:ValuesInvalidSize',nbins));
        end
        right = strcmp(validatestring(...
            p2.Results.IncludedEdge,{'left','right'}),'right');
    else
        values_provided = false;
        right = false;
    end
    
    bins = discretizemex(x, edges, right);
    if values_provided
        values(end+1) = NaN;
        bins(isnan(bins)) = length(values);
        % reshape needed when x and values are vectors of different orientation 
        bins = reshape(values(bins),size(x)); 
    end
    
end

end

function names = gencatnames(edges,right)

if right
    leftedge = '(';
    rightedge = ']';
else
    leftedge = '[';
    rightedge = ')';    
end

nbins = length(edges)-1;
names = cell(1,nbins);
for i = 1:nbins
    names{i} = sprintf([leftedge '%0.5g, %0.5g' rightedge],edges(i),edges(i+1));
end

if right
    names{1}(1) = '[';
else
    names{end}(end) = ']';
end

end


       
