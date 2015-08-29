%dati iniziali esatti

rho =  1.0e+03 *[  0.504129720000000   0.465341710184454   0.406940199998636;
  -0.514027600000000  -0.249446841939143   0.015785584115931;
   1.284749170000000   1.552664856957740   1.816111540852702];


%posizioni dell'osservatore
R = 1.0e+03 *[ 6.378137000000000   6.378094862166766   6.377969375321201; 0   0.023184476937970   0.046241088156316; 0 0  0];

r = 1.0e+03 *[ 6.882266720000000   6.843436546783947   6.784908546954059;
  -0.514027600000000  -0.226255332460717   0.062168296452558;
   1.284749170000000   1.552664856957740   1.816111540852702];

rho_hat =[ 0.342308255303443   0.283748384912019   0.218642566679389;
  -0.349028997802026  -0.152103576731988   0.008481345975778;
   0.872355327286093   0.946758340869071   0.975768156278062];


%moduli dei rho:

rho_moduli = 1.0e+03 *[1.472736085646280   1.639980119459502   1.861212142626194];

deltat =50; 
maxit=100;
sigma2=2.4240684055477e-05;
%precisione di 0.01
%l'ultimo valore comprende TUTTE le perturbazioni superiori 
contatore=zeros(1,21);
contatore2=zeros(1,21);
contatore3=zeros(1,41);
contatore4=zeros(1,41);
%si cercano le coordinate sferiche theta e phi dei rho
rho_hat2=cartesiantopolar(rho_hat);

 for i=1:2000
   
   pert=sigma2*randn(2,3);
   pert(3,:)=0; % non perturbo i moduli
   %matrice di perturbazione 
   rho_hat3=rho_hat2+pert;
   rho_hat4=polartocartesian(rho_hat3);
   [rho_est, r_est, err_jn(i), psi]=Jn(rho_moduli, R, rho_hat4, r(:,2), maxit, deltat);
   [rho_est2,r_est2,err_dc(i)]=dc(rho_hat4, R, r(:,2), deltat, maxit);
   % calcolo l'altro errore
   a=sqrt(r_est(:,2)'*r_est(:,2));
   a2=sqrt(r_est2(:,2)'*r_est2(:,2));
   errore(i)=100*( a- sqrt(r(:,2)'*r(:,2)) )/(sqrt(r(:,2)'*r(:,2)));
   errore2(i)=100*(a2-sqrt(r(:,2)'*r(:,2)))/(sqrt(r(:,2)'*r(:,2)));
   k=round(errore(i)/0.01);
   k2=round(errore2(i)/0.01);
   k3=round(err_jn(i)/0.01);
   k4=round(err_dc(i)/0.01);
      if (abs(k)<=10)
        contatore(10+k+1)=contatore(10+k+1)+1;
   
     end
     if (abs(k2)<=10)
        contatore2(10+k+1)=contatore2(10+k+1)+1;
   
     end

   
     contatore3(k3+1)=contatore3(k3+1)+1;
     contatore4(k4+1)=contatore4(k4+1)+1;
   
   
end


