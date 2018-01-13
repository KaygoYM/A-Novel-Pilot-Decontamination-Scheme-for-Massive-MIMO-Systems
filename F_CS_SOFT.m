function [P_output eu] = F_CS_SOFT(L, K, S, Beta,lambda)

eu=2;
sigma=0;
P = zeros(L, K);
beta_2=zeros(K,L);
for l=1:L
    for k=1:K
        beta_2(k,l)=Beta(k,l,l)^2;
    end
end

% for k=1:K
%     sigma=sigma+beta_2(k,1);
% end
% 
% sigma=sigma/k*lambda;
% 
% for k=1:K
%   if beta_2(k,1)<sigma
%       eu=eu+1;
%   end
% end



% for l=1:L
%     [v pos]=sort(beta_2(:,l));
%     max_pos=pos(K-eu+1:K);
%     for i=1:eu
%      P(l,max_pos(i))=-1;
%     end
% end

for l=1:L
    [v pos]=sort(beta_2(:,l));
    min_pos=pos(1:eu);
    for i=1:eu
      P(l,min_pos(i))=100+i;
    end 
end

        for i = 1:L                                                      % random pilot assignment
            temp = randperm(S);
           eum=0;
            for k=1:K
               if P(i,k) == -1
                   eum=eum+1;
               end
               if P(i,k) == 0
                    P(i,k) = temp(k-eum);
               end
            end
           
        end
        
       P_output = P;
end

