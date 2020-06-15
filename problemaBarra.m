##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO MATRICIAL (ELEMENTO DE BARRA)
## /_____/_/  |_/___/___/___/     |
##                                |    problemaBarra: script principal
##---------CICLO LECTIVO 2020----------------------------------------------------


## CONFIGURACION

clear

pkg load symbolic; # Carga de paquete que me permite hacer operaciones algebraicas

warning ('off','OctSymPy:sym:rationalapprox');

## DECLARACIONES

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

k= @(E,A,L,x) (E*A/L)*[cos(x)^2 cos(x)*sin(x) -cos(x)^2 -sin(x)*cos(x);cos(x)*sin(x) sin(x)^2 -sin(x)*cos(x) -sin(x)*cos(x);-cos(x)^2 -cos(x)*sin(x) cos(x)^2 sin(x)*cos(x);-cos(x)*sin(x) -cos(x)*sin(x) sin(x)*cos(x) sin(x)^2];

## CARACTERISTICAS DEL MATERIAL

E=1;

A=1;

## DATOS GEOMETRICOS

NODO=[0 0;0 1;0 2]; % [Xi Yi] en fila i define nodo i

ELEMENTO=[1 2 E A;2 3 E A]; % [NodoInicial NodoFinal E A] en fila i define ubicacion de la barra i y sus propiedades 

CCx=[1 0;2 0;3 0]; % [Nodo Ux] define condicion de contorno en nodo i

CCy=[1 0]; % [Nodo Uy] define condicion de contorno en nodo i 

CARGAx=[1 0;2 0;3 0]; % [Nodo Px] define carga en nodo i

CARGAy=[3 -1000]; % [Nodo Py] define carga en nodo i

## SCRIPT

for u=1:size(ELEMENTO,1)

				% CALCULOS GEOMETRICOS POR BARRA

  nodoI=ELEMENTO(u,1);

  nodoF=ELEMENTO(u,2);

  Ufx=NODO(nodoF,1);
  
  Ufy=NODO(nodoF,2);

  Uix=NODO(nodoI,1);
  
  Uiy=NODO(nodoI,2);

  ladoY=Ufy-Uiy;

  ladoX=Ufx-Uix;

  
  Largo=sqrt(ladoY^2+ladoX^2);

  keyboard
  
  try

    beta=atan(ladoY/ladoX);  # CAMBIARLO EN FUNCION DEL SENO O COSENO PARA QUE NO SE QUEJE EN BETA 90 GRADOS

  catch

    beta=0;

  end_try_catch

 

				% MATRIZ LOCAL A GLOBAL
  
  
  conectividad=[2*nodoI-1 2*nodoI 2*nodoF-1 2*nodoF];
  
  K=k(ELEMENTO(u,3),ELEMENTO(u,4),Largo,beta); % Matriz local en cuatro grados de libertad

  
  for i=1:4

    for j=1:4

      try

	KG(conectividad(i),conectividad(j))=KG(conectividad(i),conectividad(j))+K(i,j);

      catch

	KG(conectividad(i),conectividad(j))=K(i,j);
      
      end_try_catch

    endfor
  endfor

endfor

		  
				% VECTOR DE CARGAS


GL=size(KG,1);

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


				% CONDENSACION ESTATICA DE GUYAN



ccT=size(CARGAx,1)+size(CARGAy,1);

for i=1:ccT

  if isnan(P(i))==1

    j=ccT+1;
    
    while (j<=size(P,1))

      if isnan(P(j))==0
	
	KG=inversion(KG,i,j);

	P=inversion(P,i,j);

	U=inversion(U,i,j);

	try
	  registroGuyan=[registroGuyan;i j];
	catch
	  registroGuyan=[i j];
	end_try_catch

	j=size(P,1)+1;
	
      endif

      j++;
      
    endwhile
  endif
  
endfor


				% RESOLUCION




K11=KG(1:ccT,1:ccT);

K12=KG(1:ccT,ccT+1:GL);

K21=KG(ccT+1:GL,1:ccT);

K22=KG(ccT+1:GL,ccT+1:GL);

keyboard

PII=P(1:ccT,:); % Cargas conocidas

UI=U(ccT+1:GL,:); % Desplazamientos conocidos

determinante=det(K11);

if determinante==0

  UII=zeros(size(PII));

else
  
  UII=inv(K11)*(PII-K12*UI);  % Desplazamientos desconocidos

endif


PI=K21*(UII)+K22*UI; % Cargas desconocidas


				% REARMADO DE RESULTADOS

j=1;h=j;i=j;

while (j<=size(P,1))

  
  
  if isnan(P(j))==1

    P(j)=PI(i);

    i++
    
  endif

  if isnan(U(j))==1

    U(j)=UII(h);

    h++
    
  endif

  j++;
  
endwhile



for u=1:size(registroGuyan,1)

  i=registroGuyan(u,1);

  j=registroGuyan(u,2);   
  
  KG=inversion(KG,i,j);
  
  P=inversion(P,i,j);

  U=inversion(U,i,j);

endfor


				% VERIFICACION

FX=sum(P(1:2:end)) % sumatoria en X

FY=sum(P(2:2:end)) % sumatoria en Y
