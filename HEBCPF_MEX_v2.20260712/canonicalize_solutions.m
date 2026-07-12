function [Z,VBook] = canonicalize_solutions(Z,VBook)
%CANONICALIZE_SOLUTIONS Normalize global signs and deterministic ordering.

for k=1:size(Z,2)
    first_nonzero=find(abs(Z(:,k))>1e-12,1);
    if ~isempty(first_nonzero) && Z(first_nonzero,k)<0
        Z(:,k)=-Z(:,k);
    end
end
[~,order]=sortrows(round(Z'*1e12)/1e12);
Z=Z(:,order);
VBook=VBook(order,:);
end
