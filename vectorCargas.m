##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO MATRICIAL (ELEMENTO DE BARRA)
## /_____/_/  |_/___/___/___/     |
##                                |    vectorCargas: Conformado del Vector de cargas y CC
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function [P,U]=vectorCargas(GL,CARGAx,CARGAy,CCx,CCy)
  
  
  
  P=NaN(GL,1);

  for i=1:size(CARGAx,1)

    GLx=2*CARGAx(i,1)-1;

    P(GLx)=CARGAx(i,2);
    
  endfor

  for i=1:size(CARGAy,1)

    GLy=2*CARGAy(i,1);

    P(GLy)=CARGAy(i,2);

  endfor



  
				% VECTOR DE DESPLAZAMIENTOS
  U=NaN(GL,1);

  for i=1:size(CCx,1)

    GLx=2*CCx(i,1)-1;

    U(GLx)=CCx(i,2);

  endfor

  for i=1:size(CCy,1)

    GLy=2*CCy(i,1);

    U(GLy)=CCy(i,2);

  endfor

endfunction

  
