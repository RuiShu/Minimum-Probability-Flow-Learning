function show(m)
% display(m)
%
% Display double array
%
% Place in a directory called @double
% The @double directory must be in a
% directory that is on the Matlab PATH

if any(m)
  if any(inputname(1))
    fprintf('%s =\n', inputname(1));
  end
  f=9.2;
  [I,J]=size(m);
  if (I>100 | J > 30)
    fprintf('Double array (%dx%d) only partial display...\n', I, ...
            J);
    m = m(:,1:10);
    [I,J]=size(m);
  end
  if all(abs(m)<100 & (abs(m)>1e-2 | m==0))
    fprintf([repmat(['%' num2str(f) 'f'],1,J) '\n'],m');
  else
    fprintf([repmat(['%' num2str(f) 'g'],1,J) '\n'],m');
  end
end