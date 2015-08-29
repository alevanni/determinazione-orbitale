function [rho_est, r_est, err, psi]=Jn(rho0, R, rho_hat, r_vero, maxit, deltat)
 %Jn con tre osservazioni
 %R sono le posizioni dell'osservatore messe per colonne
 %rho_hat sono i versori della posizione del corpo rispetto all'osservatore
 
 %INIZIALIZZAZIONE
%suppongo di conoscere il vero punto intermedio
%r_vero=[6.842871940239587; -0.226221921024501; 1.552495618438757]*1.0e+03; %e' r2, serve per calcolare l'errore
k_vero=sqrt(r_vero'*r_vero);


 rho_est=rho0'; %si inizializza rho 
 psi=zeros(3,1); %residui
 J=zeros(3,3); %matrice jacobiana
 mu=398600.44;
 %si inizializza r
 
   for i=1:3
       r_est(:,i)=R(:, i)+rho_est(i)*rho_hat(:,i);
   end
 
 %ITERAZIONE
  for j=1:maxit
  
    %si calcolano ck e dk
       
     for i=1:3
        a=sqrt(r_est(:,i)'*r_est(:,i));
        c(i)=1/2 + (mu*deltat^2)/(4*a^3) ;
        d(i)=c(i); %li tengo entrambi per maggiore chiarezza
     end
    
    
    psi=c(2)*(R(:,1)+rho_est(1)*rho_hat(:,1)) + d(2)*(R(:,3)+rho_est(3)*rho_hat(:,3)) - (R(:,2)+rho_est(2)*rho_hat(:,2));
    %ho calcolato il residuo
    
    %ora va calcolato lo jacobiano J
      z=sqrt(r_est(:,2)'*r_est(:,2));
      cp=-(3*mu)/(4* z^5) *(rho_hat(:,2)'*r_est(:,2))*deltat^2 ;
      J(:,1)=c(2)*rho_hat(:,1);
      J(:,2)=cp*(r_est(:,1)+r_est(:,3))-rho_hat(:,2);
      J(:,3)=c(2)*rho_hat(:,3);
      
    %ora si calcola l'incremento di rho
    rho_est= rho_est - (J^(-1))*psi;
    %adesso aggiorno r
    for i=1:3
       r_est(:,i)=R(:, i)+rho_est(i)*rho_hat(:,i);
    end
    b=r_vero-r_est(:,2);
    %ora si calcola l'errore
    c=sqrt(b'*b)*100;
    err=c/k_vero;
    %e anche l'altro errore
    %commenta se non serve
    %errore=100*( sqrt(r_est(:,2)'*r_est(:,2))- k_vero)/k_vero;
    if (err<7e-2)
       j
       %break
    end 
   end
end
