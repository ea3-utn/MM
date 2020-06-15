##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO MATRICIAL (ELEMENTO DE BARRA)
## /_____/_/  |_/___/___/___/     |
##                                |    condensacionGuyan: Cond. Estatica por el Metodo de Guyan
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function [KG,P,U,fq]=condensacionGuyan(KG,P,U,fq,registroGuyan,Ni,Nf)


  
  filas=size(registroGuyan,1)+1;
  
  for u=1:filas-1

    coord=filas-u;

    i=registroGuyan(coord,Ni);

    j=registroGuyan(coord,Nf);   
    
    KG=inversion(KG,i,j);
    
    P=inversion(P,i,j);

    U=inversion(U,i,j);

    fq=inversion(fq,i,j);

  endfor
  
endfunction

  
