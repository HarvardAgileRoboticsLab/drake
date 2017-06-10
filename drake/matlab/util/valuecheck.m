function [tf,errstr] = valuecheck(val,desired_val,tol)
% [success, message] = VALUECHECK(val, desired_val, tol)
%
%   Returns a boolean and a message describing which values are unequal. If
%   no output is specified valuecheck will succeed silently or throw an
%   error.

if (nargin<3 || isempty(tol))
    tol=1e-8;
end

tf = true;
errstr = '';

% broadcast scalar expected value to correct size
if isscalar(desired_val) && ~isempty(val)
  desired_val = repmat(desired_val,size(val));
end

% check sizes
if ((length(size(val))~=length(size(desired_val))) || any(size(val)~=size(desired_val)))
  errstr = sprintf('Wrong size.  Expected %s but got %s.', mat2str(size(desired_val)), mat2str(size(val)));
  tf = false;
  if (nargout>0)
    return;
  end
end

% check NAN positions
nan_idx_val = isnan(val(:));
nan_idx_des = isnan(desired_val(:));
tf = tf & isequal(nan_idx_val,nan_idx_des);

% check numeric values
val_idx = ~nan_idx_val & ~nan_idx_des;
if any(abs(val(val_idx)-desired_val(val_idx)) > tol)
  if (ndims(val)<=2 && length(val)<=6) % for small 2-d arrays print the whole matrix
    % clean before printing
    desired_val(abs(desired_val)<tol/2)=0;
    val(abs(val)<tol/2)=0;
    errstr = sprintf('Values don''t match.  Expected \n%s\n but got \n%s', mat2str(desired_val), mat2str(val));
  else % for large matrices or N-d arrays print only the mismatched values 
    err = desired_val-val;
    ind=find(abs(err(:))>tol | (nan_idx_val(:) ~= nan_idx_des(:)));
    val_size = size(desired_val);
    a = cell(1,length(val_size));
    [a{:}] = ind2sub(val_size,ind);
    errstr_cell_array = cell(numel(ind),1);
    for i=1:numel(ind)
      b = cellfun(@(b) {num2str(b(i))},a);
      indstr = ['(',strjoin(b,','),')'];
      errstr_cell_array{i} = sprintf('%10s %12f %15f\n',indstr,full(val(ind(i))),full(desired_val(ind(i))));
    end
    errstr = sprintf('Values don''t match.\n    Index       Value       Desired Value\n   -------     --------    ---------------\n%s', [errstr_cell_array{:}]);
  end
  
  tf = false;
end

% if the values are not equal and there are no outputs throw an error
if ~tf && (nargout == 0)
  error('Drake:ValueCheck',errstr); 
end
