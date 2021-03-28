function OA = OA_CM(X, AM)

[~,in]=max(AM,[],1);
[~,in1]=max(X,[],1);
 OA = sum(in==in1)/length(in);
end
