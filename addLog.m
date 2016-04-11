function v = addLog(x)
  if length(x) == 0,
    v = -Inf;
  else
    [N,D] = size(x);
    v = log(sum(exp(x - repmat(max(x),N,1)))) + max(x);
  end;
  
%    v     = x(1,:);
%    for d=1:D,
%      for i=2:N,
%        v(d) = addLog2(v(d),x(i,d));
%      end;
%    end;
%  end;
 
 
 
% function v = addLog2(x,y)
%  zero = -Inf;
%  if x == zero, v = y;
%  elseif y == zero, v = x;
%  elseif x - y > 32, v = x;
%  elseif x > y, v = x + log(1+exp(y-x));
%  elseif y - x > 32, v = y;
%  else v = y + log(1+exp(x-y));
%  end;
