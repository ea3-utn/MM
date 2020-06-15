##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO MATRICIAL (ELEMENTO DE BARRA)
## /_____/_/  |_/___/___/___/     |
##                                |    postProcesado: Visualizacion
##                                |                     
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function postProcesado(NODO,U,FS)

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];

  for u=1:size(NODO,1)
    
    

    try
    
    DEF=[DEF;NODO(u,1)+U(u*2-1)*FS NODO(u,2)+U(u*2)*FS];

    catch

    DEF=[NODO(u,1)+U(u*2-1)*FS NODO(u,2)+U(u*2)*FS];

    end_try_catch
    
  endfor

  figure (1);clf;hold on;grid on;

  title ('METODO MATRICIAL')

  keyboard
  
  plot(NODO(:,1),NODO(:,2),["--" markStyle(1) color(1) ";MATERIAL;"]);

  plot(DEF(:,1),DEF(:,2),["--" markStyle(2) color(2) ";DEFORMADA;"]);

  hold off

    
endfunction

  
