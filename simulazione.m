format long
%posizioni dell'osservatore (in colonna)
R = 1.0e+03 *[ 6.378137000000000   6.378094862166766   6.377969375321201; 0   0.023184476937970   0.046241088156316; 0 0  0];




%posizione e velocita' iniziale da inserire in ode45
r0=[6.88226672e+03  -5.7870e-01  -5.1402760e+02   5.7434e+00   1.28474917e+03   5.3979e+00]';

%perturbazione fissata a 
sigmas=2.4240684055477e-05;
%va fatto variare tspan 
cont=1;
for deltats=20:5:120
    deltats;
    [t,sol1]=ode45(@tb, [0; deltats], r0 );
    [t,sol2]=ode45(@tb, [0; 2*deltats], r0 );
    %il sistema e' in dimensione 6, si staccano gli elementi che ci interessano
    r_trues(:,1)=[r0(1); r0(3); r0(5)]; 
    r_trues(:,2)=[sol1(size(sol1,1),1); sol1(size(sol1,1),3); sol1(size(sol1,1),5)];
    r_trues(:,3)=[sol2(size(sol2,1),1); sol2(size(sol2,1),3); sol2(size(sol2,1),5)];
    %adesso ho le posizioni iniziali esatte
    %si cercano i rho esatti (vettori)
    rho_trues=r_trues-R; 
    
    for i=1:3
      %ora se ne cercano i moduli
      rho_ms(i)=sqrt(rho_trues(:,i)'*rho_trues(:,i));
      rho_hats(:,i)=rho_trues(:,i)/rho_ms(i);
    end
    
    
    %si cercano le coordinate sferiche theta e phi dei rho
      rho_hat2s=cartesiantopolar(rho_hats);

    %si perturbano le coordinate sferiche con sigma fisso a 5''
      %matrice di perturbazione 
  
     pert=sigmas*randn(2,3);
     pert(3,:)=0; % non perturbo i moduli
     rho_hat3s=rho_hat2s+pert;
     %torno in coordinate cartesiane
     rho_hat4s=polartocartesian(rho_hat3s); %si trovano i rho_hat perturbati

     [rho_ests, r_ests, err_jns(cont), psis]=Jn(rho_ms, R, rho_hat4s, r_trues(:,2), maxit, deltats);
     [rho_est2s,r_est2s,err_dcs(cont)]=dc(rho_hat4s, R, r_trues(:,2), deltats, maxit);


     cont=cont+1;


    
end
