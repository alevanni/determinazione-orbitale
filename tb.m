function eta=tb(t,u)
%u ha seicomponenti
%l'equazione e' del secondo ordine
%ATTENZIONE, U E' FATTO COSI'!! u = [ x(t), x'(t), y(t), y'(t) , z(t), z'(t)]
%occhio all'ordine quando si passano i valori iniziali
r = sqrt ( u(1)^2 + u(3)^2 +u(5)^2);
C = 398600.4418 ;
eta(1) =   u(2) ;
eta(2) = - (C*u(1)) / r^3 ;
eta(3) =   u(4) ;
eta(4) = - (C*u(3)) / r^3 ;
eta(5) =   u(6) ;
eta(6) = - (C*u(5)) / r^3 ;
eta=eta' ;

end
%in r0 ci sono anche le velocita'
%r0=[6.88226672e+03  -5.7870e-01  -5.140276e+02   5.7434e+00   1.28474917e+03   5.3979e+00]'
%r1=r2 = 1.0e+03 *[ 6.842838558468378  -0.226219515413950   1.552544682833838];      50s
%r3 =  1.0e+03 *[  6.782163650268242   0.062287103801379   1.815515868698273]; 100s

