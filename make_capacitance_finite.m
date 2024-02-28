function C = make_capacitance_finite(N,lij)
% MAKE_CAPACITANCE construct the capacitance matrix
%   N:      number of resonators
%   lij:    spacing between the resonators

    if N == 2
        D1 = cat(2,[1/lij(1)],[1/lij(1)]);
    elseif N > 2
        D1 = cat(2,[1/lij(1)],1./lij(1:end-1)+1./lij(2:end),[1/lij(end-1)]);
    end
    D2 = -1./lij;
    C = diag(D1) + diag(D2,1) + diag(D2,-1);

end