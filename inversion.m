##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
## _______     ___       __   __   __  |
##|   ____|   /   \     |  | |  | |  | | CATEDRA ESTRUCTURAS AERONAUTICAS III
##|  |__     /  ^  \    |  | |  | |  | |
##|   __|   /  /_\  \   |  | |  | |  | | METODO MATRICIAL PROBLEMA DE BARRAS
##|  |____ /  _____  \  |  | |  | |  | |
##|_______/__/     \__\ |__| |__| |__| | inversion: inversor de grados libert.
##                                     |
##---------CICLO LECTIVO 2019----------------------------------------------------

function [I]=inversion(z,a,b)


  
  columnas=size(z,2);
  
  for i=1:size(z,1)

    if (columnas>1) % Si es matriz
      
      for j=1:columnas

	if (i==a)
	  c=b;
	elseif (i==b)
	  c=a;
	else
	  c=i;
	endif
	
	if (j==b)
	  d=a;
	elseif (j==a)
	  d=b;
	else
	  d=j;
	endif
	
	I(c,d)=z(i,j);

      endfor

    else   % Si es vector

      if (i==a)
	c=b;
      elseif (i==b)
	c=a;
      else
	c=i;
      endif

      I(c,1)=z(i,1);
      
    endif
    
      
  endfor

endfunction

