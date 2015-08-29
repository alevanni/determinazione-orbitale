
sigmap=2.4240684055477e-05;
for i=1:2000
   q=2000*rand(1,3); %perturbo il dato rho_iniziale con probabilita' uniforme

   rho_hat2=cartesiantopolar(rho_hat);
   %matrice di perturbazione 
  
    pert=sigmap*randn(2,3);
    pert(3,:)=0; % non perturbo i moduli

    rho_hat3=rho_hat2+pert;
    %torno in coordinate cartesiane
    rho_hat4=polartocartesian(rho_hat3);
    [rho_estp, r_estp, errp(i), errprel(i), psip,j(i)]=Jn(q, R, rho_hat4, r_vero, maxit, deltat);
    if (errp(i)<7e-2)
      plot3(q(1),q(2),q(3), 'o'); %se converge lo stampo
      hold on
    
    else 
        if (errp(i)<0.4)
           plot3(q(1),q(2),q(3), 'og')
        else
           plot3(q(1),q(2),q(3), 'or')
        end
    end
end
