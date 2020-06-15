##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO MATRICIAL (ELEMENTO DE BARRA)
## /_____/_/  |_/___/___/___/     |
##                                |    postProcesado: Visualizacion
##                                |                     
##---------CICLO LECTIVO 2020----------------------------------------------------------------

function postProcesado(ELEMENTO,NODO,U,FS)

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];


  for u=1:size(ELEMENTO,1)

    nodoI=ELEMENTO(u,1);

    nodoF=ELEMENTO(u,2);

    Ufx=NODO(nodoF,1);
    
    Ufy=NODO(nodoF,2);

    Uix=NODO(nodoI,1);
    
    Uiy=NODO(nodoI,2);


    try
      
      DEF=[DEF;Uix+U(nodoI*2-1)*FS Uiy+U(nodoI*2)*FS;Ufx+U(nodoF*2-1)*FS Ufy+U(nodoF*2)*FS];

      NODOplot=[NODOplot;Uix Uiy;Ufx Ufy];

    catch

      DEF=[Uix+U(nodoI*2-1)*FS Uiy+U(nodoI*2)*FS;Ufx+U(nodoF*2-1)*FS Ufy+U(nodoF*2)*FS];

      NODOplot=[Uix Uiy;Ufx Ufy];

    end_try_catch

    
  endfor

  figure (1);clf;hold on;grid on;

  title ('METODO MATRICIAL')
  
  plot(NODOplot(:,1),NODOplot(:,2),["--" markStyle(2) color(2) ";MATERIAL;"]);

  plot(DEF(:,1),DEF(:,2),["--" markStyle(4) color(4) ";DEFORMADA;"]);

  hold off

    
endfunction

  
