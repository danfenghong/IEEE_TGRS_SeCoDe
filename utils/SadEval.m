function [SAD, SADerr, perm_mtx] = SadEval(A_est, A_grt, c)
%% permute results
CRD = corrcoef([A_grt A_est]);
DD = abs(CRD(c+1:2*c,1:c));  
perm_mtx = zeros(c,c);
aux = zeros(c,1);
for i=1:c
    [ld cd] = find(max(DD(:))==DD); 
    ld=ld(1); cd=cd(1); % in the case of more than one maximum
    perm_mtx(ld,cd)=1; 
    DD(:,cd)=aux; DD(ld,:)=aux';
end
A_est = A_est*perm_mtx;

%% quantitative evaluation of spectral signature and abundance
SADerr = zeros(1,c);
for i =1:c
    SADerr(i) = acos(A_grt(:,i)'*A_est(:,i)/norm(A_grt(:,i))/norm(A_est(:,i)));
end
SAD = mean(SADerr);
end