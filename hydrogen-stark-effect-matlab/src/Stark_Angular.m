function Iang = Stark_Angular(li,lj)

% ---------------------------------------------------------
% Angular matrix element
%
% <li,0 | cos(theta) | lj,0>
%
% Nonzero only if Delta l = +-1
% ---------------------------------------------------------

if li == lj + 1

    Iang = (lj+1) / ...
        sqrt((2*lj+1)*(2*lj+3));

elseif li == lj - 1

    Iang = lj / ...
        sqrt((2*lj-1)*(2*lj+1));

else

    Iang = 0;

end

end