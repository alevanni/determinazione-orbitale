function x=gradiente_coniugato(A,b)
% Calcolo il residuo r = b -  Ax
 % che per x=0 coincide con b
  x=zeros(length(b),1);
  r=b;
  for i=1:length(b)
     rho = r'*r; %norma quadra del resto
     if (i==1) 
        p = r;
     else
        beta = rho / rhop ;
        p = r + beta * p;
     end

     
     % Calcolo q = Ap
     q=A*p;

        %  q=q+p*0.001d0  regolarizzazione: scommentare per regolarizzare

     alpha = rho/(p'*q); %prodottoscalare

     %Ci prepariamo a passare alla prossima iterazione
     x = x + alpha*p;

     % Calcolo il prossimo residuo; r = b - Ax ; ora x = x + alpha * p
     % e quindi r diventa b - A(x + alpha*p) = b - Ax - alpha*q =
     % r - alpha * q
     r = r - alpha * q;
     rhop = rho;
      
     
     
     
  end

end
