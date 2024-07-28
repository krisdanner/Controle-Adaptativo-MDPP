%--------------------------------------------
clear;
close all;
clc;
echo on

% ---------- Adaptativo MDPP (Minimum-Degree Pole Placement)----------
% Objetivo - Práticas em Controle Discreto orientadas para 
%          - Controle Digital de Nivel de Tanques 
%-----------------------------------------------------------
%  Autor: Christian Danner Ramos de Carvalho
%-----------------------------------------------------------
%    Controle de Nível de Tanques
%-----------------------------------------------------------

echo off

%---------------- Bloco 1 - Setup da Planta -------------------------
%---------------- Modelo do Tanque     ------------------------------
%---------------- Parâmetros do Tanque ---------------------------------
R1=0.5;R2=0.4;C1= 144;C2= 144;
%----------- Periodo de Amostragem e tempo total de simulação
Ts=1;t=700;
%------------Descrições:  a)ESPAÇO DE ESTADO  ---------------------------
A=[-(1/(R1*C1)),0;1/(R1*C2),-1/(R2*C2)];B=[1/C1;0]; C=[0 1];D=[0];
sys=ss(A,B,C,D);
%---                    b) FUNÇÃO DE TRANSFERÊNCIA--------------------
[num_c,den_c]=ss2tf(A,B,C,D);
G=tf(num_c,den_c);  % função de tranferencia do sistema
%------------------------------------------------------------------------
% ------------------ DISCRETIZAÇÃO da Planta (atuador+Processo)
sysd=c2d(sys,Ts,'zoh');        %Espaço de Estados
H=c2d(G,Ts,'zoh');             %Função de Transferência
[num_h,den_h] = tfdata(H);
den_h=den_h{:};num_h=num_h{:};  % Armazenas os coeficientes den e num
%-----------------------------------------------------------------------
%-------------------- Bloco 2 -   Especificações de Projeto
%--------------------------------------------------------------
% ----------------- Sistema em Malha Fechada desejado ------------------- 
zeta = 0.70;                    % Coeficiente de amortecimento
wn   = 0.0369;                  % Frequência natural não amortecida
pd   = exp(-zeta*wn*Ts +...
       j*wn*sqrt(1-zeta^2)*Ts); % Pólos discretos desejados
Am=poly([pd conj(pd)]);         % Am da FT desejada em MF
Bm=[polyval(Am,1) 0];           % Bm da FT desejada em MF
%-----------------------------------------------------------------------
AmBm_param_velho_A = [Am(2) Am(3) Bm(1)];

%------------------------------------------------------------------
%----------- Bloco 3 - Controlador Adaptativo Indireto (MDPP)
%----------------------------------------------------------------
P=diag([1 1 1 1]);              % Matriz de covariância inicial
lambda=0.5;                     % Fator de Esquecimento
%---------------
% -------------------------- Constante Inciais -------------------------- %
teta0=[0 0.001 0.1 0.1]';       % Coeficientes iniciais
%teta0=[0 0 0.1 0]';
teta = [teta0 teta0 teta0];     % Armazena teta para os primeiros instantes
y(1:3) = zeros(1,3);            % Armazena y para os primeiros instantes
u(1:3) = zeros(1,3);            % Armazena u para os primeiros instantes
[yc,t] = gensig('square',400,t,Ts); % Cria a Referêcia(yc) e o vetor t
yc=~yc*7;                       % Referencia com valor superior de 7cm
%-------------------------------------------------------------
%------------ Bloco 4  -  Malha o controle Adaptativo - MDPP-RLS
%-------------------------------------------------------------
for k=3:length(t)
% ----------------------- Simulador ------------
% a) Sendor de medição da saída - Nivel do tanque 2
 y(k)=-den_h(2)*y(k-1) -den_h(3)*y(k-2) ...    % Saída da Planta
      +num_h(2)*u(k-1) +num_h(3)*u(k-2); 
%------------------------------Bloco 5  
% ----- Estimador RLS with exponential forgetting
 Fi(:,k-1) = [-y(k-1) -y(k-2) u(k-1) u(k-2)]'; % Atualização de Fi
 eps(k) = y(k) - Fi(:,k-1)'*teta(:,k-1);       % Erro de estimação
  x = (lambda + Fi(:,k-1)' * P * Fi(:,k-1) );  % Denominador da Eq. eps
 K = P* Fi(:,k-1) / x;                         % Ganho do estimador
 teta(:,k) = teta(:,k-1) + K*eps(k);           % Vetor de estim. de parâm.
 P = (P-K*Fi(:,k-1)'*P)/lambda;                % Matriz de covariância
  na1 = teta(1,k);                     % a_1 Estimado
  na2 = teta(2,k);                     % a_2 Estimado
  nb0 = teta(3,k);                     % b_0 Estimado
  nb1 = teta(4,k);                     % b_1 Estimado
 nab_param = [na1 na2 nb0 nb1]; 
%------ Bloco 5 Determinação do parâmetros do  Controlador MDPP
% ---------Minimum-degree pole placement (MDPP) ---------------- 
   S(k)  = (Am(3)-na2)/nb0;
   R(k)  = (nb1/nb0);
   S0(k) = (Am(2)- na1)/nb0;
   T(k)  = Bm(1)/nb0;
      AmBm_param_velho_B = [Am(2) Am(3) Bm(1)];
      SRT_param_velho = [S(k) R(k) S0(k) T(k)];
%------------- Sinal do Atuador      
u(k) = -R(k)*u(k-1) +T(k)*yc(k) -S0(k)*y(k) -S(k)*y(k-1);% Lei de Controle
     if u(k)<0   % Saturação para vazão negativa
      u(k)=0;
     end
end
%---------------------------------------------------------------------
%----- Bloco 6 - Saidas - Resultados baseado em Modelos
%--------------------------------------
ifig = 0; ifig =ifig + 1; figure(ifig)
%-------------------------------- fig 1
  subplot(2,1,1)
     stairs(t,teta(1,:),'b')
      hold on
     stairs(t,teta(2,:),'r')
      hold off
     legend('a1 estimado','a2 estimado')
     title('Coeficientes Estimados')
  subplot(2,1,2)
     stairs(t,teta(3,:),'b')
      hold on
     stairs(t,teta(4,:),'r')
      hold off
     legend('b0 estimado','b1 estimado')
%      axis([0 700 0 0.01])
%------------------------------fig 2
 ifig =ifig + 1; figure(ifig)
%--------------------------------------
  subplot(2,1,1)
    stairs(t,u,'b')
    legend('u - Saída da Bomba')
      legend({'$u$'},'FontSize',14,'Interpreter','latex')
      ylabel('Vazão (cm³/s)')
      title('Lei de Controle')
  subplot(2,1,2)
    plot(t,yc,t,y) 
    legend('yc','y - H3')
      ylabel('Altura (cm)')
      xlabel('Tempo (s)')  
      title('Saída da Planta - H3')
%------------------------------fig 3
 ifig =ifig + 1; figure(ifig)
%--------------------------------------      
      figure(3)
     step(tf(Bm,Am,Ts))
     title('Saída da Planta - H3')
% publish('Simulador_Medic_Armazenam_C_RLC_serie_p.m','pdf')     