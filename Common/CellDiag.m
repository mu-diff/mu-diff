function DiagCell = CellDiag(StringValue, varargin)
    DiagCell = cell(varargin{1:end});
    ni = size(DiagCell,1);
    nj = size(DiagCell,2);
    nk = size(DiagCell,3);
    if( nk > 1)
        error('cellDiag is for 2d use only');
    end
    mindim = min(ni,nj);
    for i = 1:mindim
       DiagCell{i,i} = StringValue; 
    end
end