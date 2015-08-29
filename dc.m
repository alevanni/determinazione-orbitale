function [rho_est,r,err]=dc(rho_hat, R, r_vero, deltat, maxit)
%maxit e' il massimo numero di iterazioni
%rho_hat sono i versori della posizione del corpo rispetto all'osservatore.
%rho_hat e' una matrice, i versori sono messi in colonna
%mi faccio restituire anche r cosi' calcolo l'errore
n=size(rho_hat,2); %numero di colonne: da' il numero delle osservazioni
%ora posso riempire la matrice M^T*M che ha dimensione 3x 3
%la chiamero' C (matrice normale)
k=3*(n-2); %dimensione della matrice C
mu=398600.44; %parametro gravitazionale della terra
 %se non si mantengono separate variabili di entrata e di uscita, diverse prove in successione hanno risultati corrotti
rho_est=zeros(3,1); %l'iterazione parte da 0  
r=R;
 %r e' la matrice che contiene i vettori r delle posizioni reali del corpo, per colonne
 %inizialmente e' uguale a r, perche' rho e' zero

%suppongo di conoscere il vero punto intermedio
%r_vero=[6.842871940239587; -0.226221921024501; 1.552495618438757]*1.0e+03; %e' r2, serve per calcolare l'errore
k_vero=sqrt(r_vero'*r_vero);


for j=1:maxit
  %si calcolano i vettori c e d  
  
  for i=1:3
     a=sqrt(r(:,i)'*r(:,i));
     c(i)=(1/2+(mu*deltat^2)/(4*a^3) );
     d(i)=c(i); %li tengo entrambi per maggiore chiarezza 
  end
  
  %riempo M:
  M(:,1)=c(2)*rho_hat(:,1);
  M(:,2)=-rho_hat(:,2);
  M(:,3)=d(2)*rho_hat(:,3); 
 %riempo C 
  C(1,1)=c(2)^2;
  C(1,2)=-c(2)*rho_hat(:,1)'*rho_hat(:,2);
  C(1,3)=c(2)*d(2)*rho_hat(:,1)'*rho_hat(:,3);
  C(2,1)=C(1,2);
  C(2,2)=1;
  C(2,3)=-d(2)*rho_hat(:,2)'*rho_hat(:,3);
  C(3,1)=C(1,3);
  C(3,2)=C(2,3);
  C(3,3)=d(2)^2;
  IC=C^(-1);
  %ora va trovato xi (ho tre osservazioni, quindi un solo residuo, cioe' xi2)
  xi=R(:,2)-c(2)*R(:,1)-d(2)*R(:,3);
  %trovo il nuovo rho
  rho_est=(IC*M')*xi;
  %adesso aggiorno r
    for i=1:3
       r(:,i)=R(:, i)+rho_est(i)*rho_hat(:,i);
    end

  %ora si calcola l'errore 
    b=r_vero-r(:,2);
    %ora si calcola l'errore
    v=sqrt(b'*b)*100;
    err=v/k_vero;
    
    if (err<6e-2)
       j;
       break
    end
  end   
end
