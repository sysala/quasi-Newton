function r_crit=a_function(kappa,alpha,beta)

% this function find critical values of the function a according to a
% prescribed value kappa. The function a depends on material parameters
% alpha and beta

r=0:0.01:10;
a=alpha+(1-alpha)*r.^8./(r.^8+beta);
a_r=8*beta*(1-alpha)*(r.^7)./((r.^8+beta).^2);
a_bar=a+r.*a_r;

j=zeros(1,1001); j(1)=1;
k_hist=zeros(1,1001);
ind=1; ind2=1;
i=1;
while true
    k_i=a_bar(i)/a(i);    
    k_hist(ind)=k_i;
    if k_i<kappa
      k=k_i;  
      while k<kappa+1
        ind=ind+1;  
        k=max(a_bar(i:ind))/min(a(i:ind));        
        if ind==1001
            break
        end
        k_hist(ind)=k;
      end % while 
      if ind==1001
        break
      else
        i=ind;
        ind2=ind2+1;
        j(ind2)=i;
      end
    else
      k=k_i;  
      while k>=kappa-1
        ind=ind+1;  
        k=max(a_bar(i:ind))/min(a(i:ind));        
        if ind==1001
            break
        end
        k_hist(ind)=k;
      end % while 
      if ind==1001
        break
      else
        i=ind;
        ind2=ind2+1;
        j(ind2)=i;
      end
    end
    
end % while true

j=j(1:ind2);
r_crit=r(j);

end