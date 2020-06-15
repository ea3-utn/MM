##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
##     _________    ____________  |    
##    / ____/   |  /  _/  _/  _/  |    CATEDRA ESTRUCTURAS AERONAUTICAS III
##   / __/ / /| |  / / / / / /    |
##  / /___/ ___ |_/ /_/ /_/ /     |    METODO MATRICIAL (ELEMENTO DE BARRA)
## /_____/_/  |_/___/___/___/     |
##                                |    metodoM: script principal
##---------CICLO LECTIVO 2020----------------------------------------------------

## CONFIGURACION

clear

pkg load symbolic; # Carga de paquete que me permite hacer operaciones algebraicas

warning ('off','OctSymPy:sym:rationalapprox');

## DECLARACIONES

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

## CARACTERISTICAS DEL MATERIAL

E1=70000;
E2=200000;
A1=200;
h=3000;
A2=350;
d1=2700;
A3=150;
d2=700;

## CARACTERISTICAS FISICO-GEOMETRICAS

NODO=[0 h;d1 h;d2 0]; % [Xi Yi] en fila i define nodo i

ELEMENTO=[1 2 E1 A1;2 3 E2 A2;3 1 E1 A3]; % [NodoInicial NodoFinal E A] en fila i define ubicacion de la barra i y sus propiedades 

## CARGAS ---> L O C A L E S

cargaLocal=[1 -135+80500 -135-80500;2 168000 -168000;3 36225 -36225]; # Nelemento Ni Nf

## CONDICIONES DE CONTORNO Y CARGAS ---> G L O B A L E S

CCx=[1 0]; % [Nodo Ux] define condicion de contorno en nodo i

CCy=[1 0;3 0]; % [Nodo Uy] define condicion de contorno en nodo i 

CARGAx=[2 7000;3 0]; % [Nodo Px] define carga en nodo i

CARGAy=[2 -9000]; % [Nodo Py] define carga en nodo i


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  SCRIPT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[KG,fq]=rigidezGlobal(ELEMENTO,NODO,cargaLocal);

GL=size(KG,1); % Cant. de grados de libertad globales

[P,U]=vectorCargas(GL,CARGAx,CARGAy,CCx,CCy);

				% Guyan

INDEX=linspace(1,GL,GL);

i=1;guyan=[1 1];

CC=sum(isnan(U));

while (isempty(guyan)<1 & i<CC+1)
  
  Null=INDEX(isnan(U));

  Zeros=find(U==0);

  if (Null(i)>Zeros(1))

    guyan=[Null(i) Zeros(1)];

    try

      registroGuyan=[registroGuyan;guyan];

    catch

      registroGuyan=[guyan];

    end_try_catch

    [KG,P,U,fq]=condensacionGuyan(KG,P,U,fq,guyan,1,2);

    i++;
    
  else

    guyan=[];
    
  endif
    
endwhile


				% RESOLUCION



K11=KG(1:CC,1:CC);

K12=KG(1:CC,CC+1:GL);

K21=KG(CC+1:GL,1:CC);

K22=KG(CC+1:GL,CC+1:GL);

PII=P(1:CC,:); % Cargas conocidas

Fl=fq(1:CC,:); % Cargas locales asociadas a P conocidas

Fp=fq(CC+1:end,:); % Cargas locales NO asociadas a P conocidas

UI=U(CC+1:GL,:); % Desplazamientos conocidos

determinante=det(K11);

if determinante==0

  UII=zeros(size(PII));

else
  
  UII=inv(K11)*(PII-Fl-K12*UI);  % Desplazamientos desconocidos

endif


PI=K21*(UII)+K22*UI+Fp; % Cargas desconocidas


				% REARMADO DE RESULTADOS

j=1;h=j;i=j;

while (j<=size(P,1))

  
  
  if isnan(P(j))==1

    P(j)=PI(i);

    i++;
    
  endif

  if isnan(U(j))==1

    U(j)=UII(h);

    h++;
    
  endif

  j++;
  
endwhile


[KG,P,U,fq]=condensacionGuyan(KG,P,U,fq,registroGuyan,2,1);

postProcesado(ELEMENTO,NODO,U,10)


				% VERIFICACION

FX=sum(P(1:2:end))-sum(fq(1:2:end)) % sumatoria en X

FY=sum(P(2:2:end))-sum(fq(2:2:end)) % sumatoria en Y
