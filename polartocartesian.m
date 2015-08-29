function rho_hat=polartocartesian(rho_hat2)

%passo la matrice con i versori messi in colonna, li trasformo in coordinate cartesiane e li rimetto in colonna

  n=size(rho_hat2,2);
  rho_hat=zeros(size(rho_hat2));
  for i=1:n
    
    [rho_hat(1,i),rho_hat(2,i),rho_hat(3,i)]=sph2cart(rho_hat2(1,i),rho_hat2(2,i),rho_hat2(3,i));

  end

end
