
function[T,phi1,phi2,phi3,y,x1,x2,x3] = ACE(syn,span,smooth_method)
y = table2array(syn(:,1));
x1 = table2array(syn(:,2));
x2 = table2array(syn(:,3));
x3 = table2array(syn(:,4));
T = y./norm(y);
phi1 = zeros;
phi2 = zeros;
phi3 = zeros;
phi_sum = phi1+phi2+phi3;
esq = mean(T-phi_sum)^2;
esq_lim = 10E-4;
delta_esq = 1;
    while (delta_esq > esq_lim)
        while (delta_esq > esq_lim)
            phi1 = smooth(x1, T-(phi2+phi3),span,  smooth_method);
            phi2 = smooth(x2, T-(phi1+phi3),span,  smooth_method);
            phi3 = smooth(x3, T -(phi1+phi2),span, smooth_method);
            phi_sum = phi1 + phi2 + phi3;
            delta_esq = abs((mean(T-phi_sum)^2) - esq);
            esq = mean(T-phi_sum)^2;
        end
        phi_sum = phi1+phi2+phi3;
        T = smooth(y, phi_sum, span, smooth_method);
        T = T./ norm(T);
        delta_esq = abs((mean(T-phi_sum)^2)-esq);
        esq = mean(T-phi_sum)^2
    end
    
   subplot(2,2,1);
   scatter(y,T);
   title('Y vs Y_T_r');
   xlabel('Y');
   ylabel('Y_T_r');
   
   subplot(2,2,2);
   scatter(x1,phi1);
   title('x1 vs x1_T_r');
   xlabel('x1');
   ylabel('x1_T_r');
   
   subplot(2,2,3);
   scatter(x2,phi2);
   title('x2 vs x2_T_r');
   xlabel('x2');
   ylabel('x2_T_r');
   
   subplot(2,2,4);
   scatter(x3,phi3);
   title('x3 vs x3_T_r');
   xlabel('x3');
   ylabel('x3_T_r');
end



