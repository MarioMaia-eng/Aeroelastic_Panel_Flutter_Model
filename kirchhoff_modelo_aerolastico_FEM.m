%% main_kirchhoff_unimorph_eletromecanico_RL.m
% Modelo FEM eletromecanico de uma placa unimorph piezoeletrica
% Elemento retangular de Kirchhoff, 4 nos, 3 GDL mecanicos por no.
% Caso de validacao: De Marqui Jr., Erturk e Inman (2009).


clear; close all; clc;

% Controle de execucao. As FRFs por aceleracao de base do modelo original
% ficam preservadas apenas para validacao eletromecanica complementar.
executar_FRFs_base_originais = false;

%% ================================================================
% 1) Geometria do caso unimorph do artigo
% ================================================================

L  = 0.1;          % comprimento [m]
B  = 0.02;         % largura [m]
hs = 0.5e-3;       % espessura subestrutura [m]
hp = 0.4e-3;       % espessura PZT [m]

%% ================================================================
% 2) Malha
% ================================================================

nelx = 50;
nely = 10;

nnx = nelx + 1;
nny = nely + 1;

nNodes = nnx*nny;
nElem  = nelx*nely;

dx = L/nelx;
dy = B/nely;

[xGrid, yGrid] = meshgrid(linspace(0,L,nnx), linspace(0,B,nny));

nodes = zeros(nNodes,2);
id = 0;
for j = 1:nny
    for i = 1:nnx
        id = id + 1;
        nodes(id,:) = [xGrid(j,i), yGrid(j,i)];
    end
end

elements = zeros(nElem,4);
e = 0;
for j = 1:nely
    for i = 1:nelx
        n1 = (j-1)*nnx + i;
        n2 = n1 + 1;
        n4 = j*nnx + i;
        n3 = n4 + 1;
        e = e + 1;
        elements(e,:) = [n1 n2 n3 n4];
    end
end

%% ================================================================
% 3) Materiais
% ================================================================

Es   = 100e9;
nus  = 0.30;
rhos = 7165;

Qs = Es/(1 - nus^2)*[1 nus 0; nus 1 0; 0 0 (1-nus)/2];

rhop = 7800;

c11 = 120.3e9;
c22 = 120.3e9;
c12 = 75.2e9;
c13 = 75.1e9;
c23 = 75.1e9;
c33 = 110.9e9;
c66 = 22.7e9;

% Constantes piezoeletricas 3D [C/m^2]
% Tabela do artigo: e31 = e32 = -5.2, e33 = 15.9.
% Se a sua convencao anterior usava +5.2, altere e31_3D/e32_3D aqui.
e31_3D = -5.2;
e32_3D = -5.2;
e33_3D = 15.9;

eps0 = 8.854187817e-12;
eps33S_3D = 1800*eps0;

% Reducao para estado plano de tensoes, Eq. (26)
c11b = c11 - c13^2/c33;
c22b = c22 - c23^2/c33;
c12b = c12 - c13*c23/c33;
c66b = c66;

Qp = [c11b c12b 0; c12b c22b 0; 0 0 c66b];

e31b = e31_3D - c13*e33_3D/c33;
e32b = e32_3D - c23*e33_3D/c33;
eps33b = eps33S_3D + e33_3D^2/c33;

fprintf('\nPropriedades piezoeletricas reduzidas:\n');
fprintf('e31b = %.6e C/m^2\n', e31b);
fprintf('e32b = %.6e C/m^2\n', e32b);
fprintf('eps33b = %.6e F/m\n', eps33b);

%% ================================================================
% 4) Camadas e linha neutra pela secao transformada
% ================================================================

layers_raw(1).name  = 'Substrato';
layers_raw(1).z_bot = 0;
layers_raw(1).z_top = hs;
layers_raw(1).Q     = Qs;
layers_raw(1).rho   = rhos;
layers_raw(1).isPiezo  = false;
layers_raw(1).e31bar   = 0;
layers_raw(1).e32bar   = 0;
layers_raw(1).eps33bar = 0;

layers_raw(2).name  = 'PZT-5A';
layers_raw(2).z_bot = hs;
layers_raw(2).z_top = hs + hp;
layers_raw(2).Q     = Qp;
layers_raw(2).rho   = rhop;
layers_raw(2).isPiezo  = true;
layers_raw(2).e31bar   = e31b;
layers_raw(2).e32bar   = e32b;
layers_raw(2).eps33bar = eps33b;

Ss_mat = inv(Qs);
Sp_mat = inv(Qp);
Ys_eff = 1/Ss_mat(1,1);
Yp_eff = 1/Sp_mat(1,1);

[zNeutral, n_mod, hpa, hsa, hpc, ha, hb, hc] = ...
    compute_neutral_axis_transformed_section(hs, hp, Ys_eff, Yp_eff);

fprintf('\nLinha neutra pela secao transformada:\n');
fprintf('Ys_eff = %.6e Pa\n', Ys_eff);
fprintf('Yp_eff = %.6e Pa\n', Yp_eff);
fprintf('n = Ys/Yp = %.6f\n', n_mod);
fprintf('zNeutral medido da face inferior = %.6e m\n', zNeutral);

layers = layers_raw;
layers(1).z_bot = ha;
layers(1).z_top = hb;
layers(2).z_bot = hb;
layers(2).z_top = hc;

fprintf('\nPosicoes z das camadas medidas a partir da linha neutra:\n');
for k = 1:length(layers)
    fprintf('%s: z_bot = %.6e m | z_top = %.6e m\n', ...
        layers(k).name, layers(k).z_bot, layers(k).z_top);
end

%% ================================================================
% 5) Montagem global: M, K, Fbase, Htilde, Cp e piston theory
% ================================================================

ndof = 3*nNodes;
K = sparse(ndof,ndof);
M = sparse(ndof,ndof);
Fbase = sparse(ndof,1);
Htilde = sparse(ndof,1);

% Matrizes geometricas da piston theory, independentes da velocidade:
%   Gx = int_A Nw^T*(dNw/dx) dA
%   Gv = int_A Nw^T*Nw dA
% O escoamento e considerado na direcao global +x.
Gx_aero = sparse(ndof,ndof);
Gv_aero = sparse(ndof,ndof);

Cp_total = 0;
ngp = 5;

for e = 1:nElem
    conn = elements(e,:);
    edofs = zeros(1,12);
    for a = 1:4
        node = conn(a);
        edofs(3*a-2:3*a) = [3*node-2, 3*node-1, 3*node];
    end

    [Ke, Me] = kirchhoff_rect_element(dx, dy, layers, ngp);
    Fe_base = kirchhoff_base_force_element(dx, dy, layers, ngp);
    [He, Cpe] = kirchhoff_piezo_element(dx, dy, layers, ngp);
    [Gxe, Gve] = kirchhoff_piston_element(dx, dy, ngp);

    K(edofs,edofs) = K(edofs,edofs) + Ke;
    M(edofs,edofs) = M(edofs,edofs) + Me;
    Fbase(edofs) = Fbase(edofs) + Fe_base;
    Htilde(edofs) = Htilde(edofs) + He;
    Gx_aero(edofs,edofs) = Gx_aero(edofs,edofs) + Gxe;
    Gv_aero(edofs,edofs) = Gv_aero(edofs,edofs) + Gve;
    Cp_total = Cp_total + Cpe;
end

fprintf('\nCapacitancia total Cp_total = %.6e F\n', Cp_total);

%% ================================================================
% 6) Condicao de contorno: engaste em x = 0
% ================================================================

tol = 1e-12;
fixedNodes = find(abs(nodes(:,1)) < tol);
fixedDofs = [];
for i = 1:length(fixedNodes)
    node = fixedNodes(i);
    fixedDofs = [fixedDofs, 3*node-2, 3*node-1, 3*node];
end

allDofs  = 1:ndof;
freeDofs = setdiff(allDofs, fixedDofs);

Kff = K(freeDofs, freeDofs);
Mff = M(freeDofs, freeDofs);
Fbase_f = Fbase(freeDofs);
Hff = Htilde(freeDofs);
Gxff = Gx_aero(freeDofs, freeDofs);
Gvff = Gv_aero(freeDofs, freeDofs);

%% ================================================================
% 6.1) Verificacao da excitacao de base
% ================================================================

m_area_check = 0;
for k = 1:length(layers)
    h_layer = layers(k).z_top - layers(k).z_bot;
    m_area_check = m_area_check + layers(k).rho*h_layer;
end

m_total_aprox = m_area_check*L*B;

fprintf('\nVerificacao da excitacao de base:\n');
fprintf('m_area = %.6e kg/m^2\n', m_area_check);
fprintf('m_total aproximada = %.6e kg\n', m_total_aprox);
s = sum(Fbase(1:3:end));    % soma (pode ser sparse)
s = full(s);                % converte apenas o resultado para full

if isscalar(s)
    fprintf('sum(Fbase nos GDL w) = %.6e\n', s);
else
    % imprime todos os elementos (cada um em nova linha)
    fprintf('sum(Fbase nos GDL w) = %.6e\n', s);
end
fprintf('norm(Fbase) = %.6e\n', norm(Fbase));
fprintf('norm(Fbase_f) = %.6e\n', norm(Fbase_f));

%% ================================================================
% 7) Modal: curto-circuito e circuito aberto
% ================================================================

nModes = 6;
[Vsc, Dsc] = eig(full(Kff), full(Mff));
lambda_sc = diag(Dsc);
valid_sc = lambda_sc > 1e-6;
lambda_sc = lambda_sc(valid_sc);
Vsc = Vsc(:,valid_sc);
[lambda_sc, idx_sc] = sort(lambda_sc);
Vsc = Vsc(:,idx_sc);
omega_sc = sqrt(lambda_sc);
freq_sc  = omega_sc/(2*pi);

fprintf('\nFrequencias naturais em curto-circuito aproximado:\n');
for i = 1:nModes
    fprintf('Modo %d: %.4f Hz\n', i, freq_sc(i));
end

Kff_open = Kff + (Hff*Hff')/Cp_total;
[Voc, Doc] = eig(full(Kff_open), full(Mff));
lambda_oc = diag(Doc);
valid_oc = lambda_oc > 1e-6;
lambda_oc = lambda_oc(valid_oc);
Voc = Voc(:,valid_oc);
[lambda_oc, idx_oc] = sort(lambda_oc);
Voc = Voc(:,idx_oc);
omega_oc = sqrt(lambda_oc);
freq_oc  = omega_oc/(2*pi);

fprintf('\nFrequencias naturais em circuito aberto aproximado:\n');
for i = 1:nModes
    fprintf('Modo %d: %.4f Hz\n', i, freq_oc(i));
end

V = Vsc;
freq = freq_sc;
omega = omega_sc;

%% ================================================================
% 8) Amortecimento proporcional
% ================================================================

alpha_damp = 4.886;          % [rad/s]
beta_damp  = 1.2433e-5;      % [s/rad]
Cff = alpha_damp*Mff + beta_damp*Kff;

fprintf('\nAmortecimento modal estimado em curto-circuito:\n');
for m = 1:min(nModes,length(freq_sc))
    wm = 2*pi*freq_sc(m);
    zeta_m = 0.5*(alpha_damp/wm + beta_damp*wm);
    fprintf('Modo %d: f = %.4f Hz | zeta = %.6f\n', m, freq_sc(m), zeta_m);
end

%% ================================================================
% 9) PISTON THEORY: ESTABILIDADE AEROELASTICA E CIRCUITO RESISTIVO
% ================================================================
% Esta e a etapa principal para localizar a faixa de flutter.
%
% Convencao aerodinamica adotada:
%
%   Delta p = -lambda_a*w_x - mu_a*w_t
%
% aplicada no lado direito da equacao mecanica. Ao transportar a carga
% aerodinamica para o lado esquerdo, obtem-se:
%
%   M*qdd + (C + Ca)*qd + (K + Ka)*q - H*v = 0
%
%   Ka = lambda_a*Gx
%   Ca = mu_a*Gv
%
% com a piston theory linear de primeira ordem:
%
%   lambda_a = 2*q_inf/sqrt(Mach^2 - 1)
%   mu_a     = lambda_a/U * (Mach^2 - 2)/(Mach^2 - 1)
%
% Para a primeira formulacao, adota-se Mach >= 2. A aceleracao de base
% nao participa da determinacao do flutter.

executar_estabilidade_piston = true;

if executar_estabilidade_piston

    %% ------------------------------------------------------------
    % 9.1) Parametros do escoamento e da reducao modal
    % -------------------------------------------------------------
    rho_inf = 1.225;             % densidade do ar [kg/m^3]
    a_inf   = 343.0;             % velocidade do som [m/s]

    % Faixa refinada solicitada para localizar a fronteira de flutter.
    Uvec    = linspace(1000.0, 2000.0, 301);  % [m/s]
    MachVec = Uvec/a_inf;
    nU      = numel(Uvec);

    % Numero de modos estruturais usados na reducao.
    nModesAero = min(24, size(Vsc,2));
    nModesPlot = min(8, nModesAero);

    % Diagnostico modal no Command Window.
    imprimir_modos_console = true;
    nModesConsole = min(8,nModesAero);
    R_referencia_modos_FRF = 1e4;  % tabela modal detalhada em cada figura de FRF

    % Resistencias finitas entre os limites exatos de curto-circuito
    % (R_L = 0) e circuito aberto (R_L -> infinito). Os dois limites
    % exatos sao calculados separadamente para evitar aproximacoes
    % numericamente mal condicionadas.
    R_values = 10.^(1:8);             % 10^1 ate 10^8 Ohm, uma curva por decada
    nR_aero  = numel(R_values);

    fprintf('\n============================================================\n');
    fprintf('PISTON THEORY - ESTABILIDADE E CIRCUITO RESISTIVO\n');
    fprintf('============================================================\n');
    fprintf('Faixa de Mach: %.3f ate %.3f\n', MachVec(1), MachVec(end));
    fprintf('Faixa de velocidade: %.3f ate %.3f m/s\n', Uvec(1), Uvec(end));
    fprintf('Modos retidos: %d\n', nModesAero);

    %% ------------------------------------------------------------
    % 9.2) Base modal de curto-circuito e projecao das matrizes
    % -------------------------------------------------------------
    PhiA = Vsc(:,1:nModesAero);

    % Normalizacao pela massa.
    for im = 1:nModesAero
        modal_mass = real(PhiA(:,im)'*Mff*PhiA(:,im));
        if modal_mass <= 0
            error('Massa modal nao positiva no modo %d.', im);
        end
        PhiA(:,im) = PhiA(:,im)/sqrt(modal_mass);
    end

    Mr  = full(PhiA'*Mff*PhiA);
    Kr  = full(PhiA'*Kff*PhiA);
    Cr  = full(PhiA'*Cff*PhiA);
    Hr  = full(PhiA'*Hff);
    Gxr = full(PhiA'*Gxff*PhiA);
    Gvr = full(PhiA'*Gvff*PhiA);

    fprintf('||Mr-I|| = %.6e\n', norm(Mr-eye(nModesAero)));
    fprintf('||Gxr||  = %.6e\n', norm(Gxr));
    fprintf('||Gvr||  = %.6e\n', norm(Gvr));
    fprintf('Assimetria relativa de Gxr = %.6e\n', ...
        norm(Gxr-Gxr','fro')/max(norm(Gxr,'fro'),eps));

    %% ------------------------------------------------------------
    % 9.3) Coeficientes aerodinamicos ao longo da velocidade
    % -------------------------------------------------------------
    lambdaA_vec = zeros(1,nU);
    muA_vec     = zeros(1,nU);

    for iu = 1:nU
        Mach = MachVec(iu);
        U    = Uvec(iu);
        betaM = sqrt(Mach^2-1);
        qinf  = 0.5*rho_inf*U^2;

        lambdaA_vec(iu) = 2*qinf/betaM;
        muA_vec(iu) = (lambdaA_vec(iu)/U)* ...
            ((Mach^2-2)/(Mach^2-1));
    end

    %% ------------------------------------------------------------
    % 9.4) Caso AE: estabilidade sem realimentacao eletrica
    % -------------------------------------------------------------
    sigmaMax_AE = NaN(1,nU);
    freqDom_AE  = NaN(1,nU);
    gDom_AE     = NaN(1,nU);

    sigmaTrack_AE = NaN(nModesPlot,nU);
    freqTrack_AE  = NaN(nModesPlot,nU);
    gTrack_AE     = NaN(nModesPlot,nU);

    prevLam = [];
    prevVec = [];

    for iu = 1:nU
        Ka = lambdaA_vec(iu)*Gxr;
        Ca = muA_vec(iu)*Gvr;

        Kt = Kr + Ka;
        Ct = Cr + Ca;

        Aae = [zeros(nModesAero), eye(nModesAero); ...
              -(Mr\Kt),          -(Mr\Ct)];

        [Vec,D] = eig(Aae);
        lam = diag(D);

        [sigmaMax_AE(iu), freqDom_AE(iu), gDom_AE(iu)] = ...
            dominant_oscillatory_mode(lam);

        [lamSel, vecSel] = track_oscillatory_modes( ...
            lam, Vec, nModesAero, nModesPlot, prevLam, prevVec);

        if ~isempty(lamSel)
            nFound = numel(lamSel);
            sigmaTrack_AE(1:nFound,iu) = real(lamSel);
            freqTrack_AE(1:nFound,iu) = abs(imag(lamSel))/(2*pi);
            gTrack_AE(1:nFound,iu) = 2*real(lamSel)./max(abs(imag(lamSel)),eps);
            prevLam = lamSel;
            prevVec = vecSel;
        end
    end

    [Ucrit_AE, MachCrit_AE, fcrit_AE, idxCrit_AE] = ...
        find_flutter_crossing(Uvec, sigmaMax_AE, freqDom_AE, a_inf);

    fprintf('\nCaso aeroelastico sem circuito:\n');
    if isnan(Ucrit_AE)
        fprintf('Nao foi encontrado cruzamento Re(lambda)=0 na faixa analisada.\n');
        fprintf('Nao altere o sinal por tentativa. Primeiro aumente a faixa de Mach,\n');
        fprintf('verifique convergencia modal e compare com um caso de literatura.\n');
    else
        fprintf('Ucrit_AE    = %.6f m/s\n', Ucrit_AE);
        fprintf('MachCrit_AE = %.6f\n', MachCrit_AE);
        fprintf('fcrit_AE    = %.6f Hz\n', fcrit_AE);
    end

    %% ------------------------------------------------------------
    % 9.4.1) Limite exato de circuito aberto
    % -------------------------------------------------------------
    % Para circuito aberto, a equacao eletrica integrada fornece
    % v = -(H_r''*q)/C_p. Ao substituir na equacao mecanica, aparece
    % a rigidez eletromecanica adicional H_r*H_r''/C_p.
    sigmaMax_OC = NaN(1,nU);
    freqDom_OC  = NaN(1,nU);
    gDom_OC     = NaN(1,nU);

    Koc_eletromec = (Hr*Hr')/Cp_total;

    for iu = 1:nU
        Ka = lambdaA_vec(iu)*Gxr;
        Ca = muA_vec(iu)*Gvr;

        Kt_oc = Kr + Ka + Koc_eletromec;
        Ct_oc = Cr + Ca;

        Aoc = [zeros(nModesAero), eye(nModesAero); ...
              -(Mr\Kt_oc),       -(Mr\Ct_oc)];

        lam_oc = eig(Aoc);
        [sigmaMax_OC(iu), freqDom_OC(iu), gDom_OC(iu)] = ...
            dominant_oscillatory_mode(lam_oc);
    end

    [Ucrit_OC, MachCrit_OC, fcrit_OC] = ...
        find_flutter_crossing(Uvec, sigmaMax_OC, freqDom_OC, a_inf);

    fprintf('\nLimite de circuito aberto:\n');
    if isnan(Ucrit_OC)
        fprintf('Nao foi encontrado flutter de circuito aberto na faixa analisada.\n');
    else
        fprintf('Ucrit_OC    = %.6f m/s\n', Ucrit_OC);
        fprintf('MachCrit_OC = %.6f\n', MachCrit_OC);
        fprintf('fcrit_OC    = %.6f Hz\n', fcrit_OC);
    end

    %% ------------------------------------------------------------
    % 9.5) Caso AEE: circuito puramente resistivo
    % -------------------------------------------------------------
    sigmaMax_AEE = NaN(nR_aero,nU);
    freqDom_AEE  = NaN(nR_aero,nU);
    gDom_AEE     = NaN(nR_aero,nU);

    Ucrit_AEE    = NaN(nR_aero,1);
    MachCrit_AEE = NaN(nR_aero,1);
    fcrit_AEE    = NaN(nR_aero,1);

    for ir = 1:nR_aero
        Rl = R_values(ir);
        fprintf('Varredura AEE: R_L = %.4e Ohm\n', Rl);

        for iu = 1:nU
            Ka = lambdaA_vec(iu)*Gxr;
            Ca = muA_vec(iu)*Gvr;

            Kt = Kr + Ka;
            Ct = Cr + Ca;

            Aaee = [zeros(nModesAero), eye(nModesAero), zeros(nModesAero,1); ...
                   -(Mr\Kt),          -(Mr\Ct),          Mr\Hr; ...
                    zeros(1,nModesAero), -(Hr'/Cp_total), -1/(Rl*Cp_total)];

            lam = eig(Aaee);
            [sigmaMax_AEE(ir,iu), freqDom_AEE(ir,iu), gDom_AEE(ir,iu)] = ...
                dominant_oscillatory_mode(lam);
        end

        [Ucrit_AEE(ir), MachCrit_AEE(ir), fcrit_AEE(ir)] = ...
            find_flutter_crossing(Uvec, sigmaMax_AEE(ir,:), ...
                                  freqDom_AEE(ir,:), a_inf);
    end

    fprintf('\nResumo do circuito resistivo:\n');
    fprintf(' R_L [Ohm]        Ucrit [m/s]      Machcrit       fcrit [Hz]\n');
    for ir = 1:nR_aero
        fprintf('%12.4e   %14.6f   %12.6f   %12.6f\n', ...
            R_values(ir), Ucrit_AEE(ir), MachCrit_AEE(ir), fcrit_AEE(ir));
    end

    %% ------------------------------------------------------------
    % 9.6) Graficos de estabilidade sem circuito
    % -------------------------------------------------------------
    figure('Name','Estabilidade aeroelastica - piston theory');

    subplot(2,1,1); hold on;
    for im = 1:nModesPlot
        plot(Uvec, gTrack_AE(im,:), 'LineWidth', 1.15);
    end
    yline(0,'k:','LineWidth',1.1);
    if ~isnan(Ucrit_AE)
        xline(Ucrit_AE,'k--',sprintf('U_{cr}=%.1f m/s',Ucrit_AE), ...
            'LabelVerticalAlignment','bottom');
    end
    grid on;
    xlim([1000 2000]);
    ylim([-0.1 0.5]);
    xlabel('Velocidade do escoamento U_\infty [m/s]');
    ylabel('g_i = 2Re(\lambda_i)/|Im(\lambda_i)|');
    title('Evolucao do amortecimento aeroelastico');

    subplot(2,1,2); hold on;
    for im = 1:nModesPlot
        plot(Uvec, freqTrack_AE(im,:), 'LineWidth', 1.15);
    end
    if ~isnan(Ucrit_AE)
        xline(Ucrit_AE,'k--');
    end
    grid on;
    xlim([1000 2000]);
    xlabel('Velocidade do escoamento U_\infty [m/s]');
    ylabel('Frequencia [Hz]');
    title('Evolucao das frequencias aeroelasticas');

    % Impressao dos modos associados ao grafico de estabilidade.
    if exist('imprimir_modos_console','var') && imprimir_modos_console
        idxConsole = [1 nU];
        if ~isnan(idxCrit_AE)
            idxConsole = [idxConsole max(1,min(nU,idxCrit_AE))];
        end
        idxConsole = unique(idxConsole,'stable');

        fprintf('\n============================================================\n');
        fprintf('MODOS DO GRAFICO DE ESTABILIDADE AEROELASTICA\n');
        fprintf('============================================================\n');
        for kk = 1:numel(idxConsole)
            iuPrint = idxConsole(kk);
            KaPrint = lambdaA_vec(iuPrint)*Gxr;
            CaPrint = muA_vec(iuPrint)*Gvr;
            Aprint = [zeros(nModesAero), eye(nModesAero); ...
                     -(Mr\(Kr+KaPrint)), -(Mr\(Cr+CaPrint))];
            labelPrint = sprintf('AE sem circuito | U = %.3f m/s | Mach = %.3f', ...
                Uvec(iuPrint),MachVec(iuPrint));
            print_state_modes(labelPrint,Aprint,nModesAero,nModesConsole);
        end
    end

    %% ------------------------------------------------------------
    % 9.7) Grafico combinado: curto-circuito ate circuito aberto
    % -------------------------------------------------------------
    % O caso AE sem realimentacao eletrica corresponde ao limite exato
    % de curto-circuito, pois a tensao nos eletrodos e nula. As curvas
    % coloridas representam resistencias finitas, e a curva tracejada-
    % pontilhada representa o limite exato de circuito aberto.

    figure('Name','Amortecimento: curto-circuito a circuito aberto');
    hold on;

    hSC = plot(Uvec,gDom_AE,'k--','LineWidth',2.0, ...
        'DisplayName','Curto-circuito (R_L = 0)');

    % Curvas para as resistencias finitas de 10^1 a 10^8 Ohm.
    % Cada curva aparece explicitamente na legenda, sem barra de cores.
    cmapR = lines(nR_aero);
    hR = gobjects(nR_aero,1);
    for ir = 1:nR_aero
        hR(ir) = plot(Uvec,gDom_AEE(ir,:), ...
            'Color',cmapR(ir,:), ...
            'LineWidth',1.25, ...
            'DisplayName',sprintf('R_L = 1e+%02d \\Omega',round(log10(R_values(ir)))));
    end

    hOC = plot(Uvec,gDom_OC,'k-.','LineWidth',2.0, ...
        'DisplayName','Circuito aberto (R_L \rightarrow \infty)');

    yline(0,'k:','LineWidth',1.2,'HandleVisibility','off');

    % Marcadores das fronteiras exatas dos dois limites eletricos.
    if ~isnan(Ucrit_AE)
        plot(Ucrit_AE,0,'ko','MarkerFaceColor','w', ...
            'MarkerSize',7,'HandleVisibility','off');
    end
    if ~isnan(Ucrit_OC)
        plot(Ucrit_OC,0,'ks','MarkerFaceColor','k', ...
            'MarkerSize',7,'HandleVisibility','off');
    end

    grid on;
    box on;
    xlim([1000 2000]);
    ylim([-0.1 0.5]);
    xlabel('Velocidade do escoamento U_\infty [m/s]');
    ylabel('g do modo mais critico');
    title('Efeito de R_L sobre o amortecimento aeroeletroelastico');
    % Legenda convencional, no mesmo estilo do exemplo enviado.
    legend([hSC; hR; hOC], ...
        'Location','southeast', ...
        'NumColumns',1, ...
        'FontSize',9, ...
        'Interpreter','tex');

    if exist('imprimir_modos_console','var') && imprimir_modos_console
        fprintf('\n============================================================\n');
        fprintf('MODO CRITICO DO GRAFICO R_L x AMORTECIMENTO\n');
        fprintf('============================================================\n');

        if isfinite(Ucrit_AE)
            [~,iuCrit] = min(abs(Uvec-Ucrit_AE));
            KaCrit = lambdaA_vec(iuCrit)*Gxr;
            CaCrit = muA_vec(iuCrit)*Gvr;
            Acrit = [zeros(nModesAero), eye(nModesAero); ...
                     -(Mr\(Kr+KaCrit)), -(Mr\(Cr+CaCrit))];
            print_dominant_state_mode('Curto-circuito (R_L = 0)',Acrit,nModesAero,Ucrit_AE);
        end

        for irPrint = 1:nR_aero
            if ~isfinite(Ucrit_AEE(irPrint))
                fprintf('R_L = %.3e Ohm | flutter nao encontrado na faixa analisada.\n',R_values(irPrint));
                continue;
            end
            [~,iuCrit] = min(abs(Uvec-Ucrit_AEE(irPrint)));
            KaCrit = lambdaA_vec(iuCrit)*Gxr;
            CaCrit = muA_vec(iuCrit)*Gvr;
            KtCrit = Kr + KaCrit;
            CtCrit = Cr + CaCrit;
            Rcrit = R_values(irPrint);
            Acrit = [zeros(nModesAero), eye(nModesAero), zeros(nModesAero,1); ...
                     -(Mr\KtCrit), -(Mr\CtCrit), Mr\Hr; ...
                     zeros(1,nModesAero), -(Hr'/Cp_total), -1/(Rcrit*Cp_total)];
            labelCrit = sprintf('R_L = %.3e Ohm',Rcrit);
            print_dominant_state_mode(labelCrit,Acrit,nModesAero,Ucrit_AEE(irPrint));
        end

        if isfinite(Ucrit_OC)
            [~,iuCrit] = min(abs(Uvec-Ucrit_OC));
            KaCrit = lambdaA_vec(iuCrit)*Gxr;
            CaCrit = muA_vec(iuCrit)*Gvr;
            Acrit = [zeros(nModesAero), eye(nModesAero); ...
                     -(Mr\(Kr+KaCrit+Koc_eletromec)), -(Mr\(Cr+CaCrit))];
            print_dominant_state_mode('Circuito aberto',Acrit,nModesAero,Ucrit_OC);
        end
    end

    %% ------------------------------------------------------------
    % 9.8) Velocidade e frequencia de flutter versus resistencia
    % -------------------------------------------------------------
    figure('Name','Flutter em funcao da resistencia');
    subplot(2,1,1);
    semilogx(R_values,Ucrit_AEE,'o-','LineWidth',1.4,'MarkerSize',6);
    hold on;
    if ~isnan(Ucrit_AE)
        yline(Ucrit_AE,'k--','AE sem circuito');
    end
    grid on;
    xlabel('Resistencia R_L [\Omega]');
    ylabel('U_{cr} [m/s]');
    title('Velocidade critica de flutter');

    subplot(2,1,2);
    semilogx(R_values,fcrit_AEE,'o-','LineWidth',1.4,'MarkerSize',6);
    hold on;
    if ~isnan(fcrit_AE)
        yline(fcrit_AE,'k--','AE sem circuito');
    end
    grid on;
    xlabel('Resistencia R_L [\Omega]');
    ylabel('f_{cr} [Hz]');
    title('Frequencia na fronteira de flutter');

    %% ------------------------------------------------------------
    % 9.9) Mapa de estabilidade
    % -------------------------------------------------------------
    figure('Name','Mapa de estabilidade U-R');
    [XR,YU] = meshgrid(log10(R_values),Uvec);
    contourf(XR,YU,sigmaMax_AEE',24,'LineColor','none');
    hold on;
    contour(XR,YU,sigmaMax_AEE',[0 0],'k','LineWidth',2.0);
    colorbar;
    grid on;
    xlabel('log_{10}(R_L [\Omega])');
    ylabel('Velocidade do escoamento U_\infty [m/s]');
    ylim([1000 2000]);
    title('Mapa de estabilidade: max Re(\lambda)');

    %% ------------------------------------------------------------
    % 9.10) FRFs PIEZOAEROELASTICAS PARA DIFERENTES R_L
    % -------------------------------------------------------------
    % As FRFs abaixo usam aceleracao harmonica de base unitaria apenas
    % como entrada de referencia. A aerodinamica altera Kt e Ct para cada
    % velocidade de escoamento. Sao calculadas:
    %
    %   |w_rel(L)/Y0|          deslocamento relativo de ponta
    %   |V/a_b|                tensao por aceleracao de base
    %   |I/a_b|                corrente por aceleracao de base
    %   P_media/a_b^2          potencia media normalizada
    %
    % Todas as magnitudes sao apresentadas em base logaritmica:
    %   amplitudes -> 20*log10(|H|)
    %   potencia   -> 10*log10(P)
    % e o eixo de frequencia tambem e logaritmico.
    %
    % IMPORTANTE: uma FRF estacionaria so possui interpretacao fisica
    % para um sistema estavel. Por isso, a configuracao padrao mascara os
    % pontos acima da fronteira de flutter de cada resistencia.

    executar_FRFs_piezoaero = true;
    mascarar_regiao_instavel_FRF = true;
    executar_envelopes_FRF_velocidade = true;

    % ---------------------------------------------------------
    % CONFIGURACAO MANUAL DAS VELOCIDADES DAS FRFs
    % ---------------------------------------------------------
    % As velocidades abaixo sao escolhidas diretamente pelo usuario e
    % nao dependem da velocidade critica calculada pelo programa.
    % Para alterar os casos analisados, modifique somente esta linha.
    U_FRF_values_manual = 400:250:1400;  % [m/s] -> 400, 650, 900, 1150, 1400

    if executar_FRFs_piezoaero

        % ---------------------------------------------------------
        % 9.10.1) Configuracao da faixa de frequencia e velocidades
        % ---------------------------------------------------------
        % Faixa das FRFs limitada a 2500 Hz.
        freqVec_AERO = logspace(log10(10),log10(2500),1800); % [Hz]
        Omega_AERO   = 2*pi*freqVec_AERO;
        nFreq_AERO   = numel(freqVec_AERO);

        % Velocidades escolhidas manualmente pelo usuario.
        U_FRF_values = unique(round(U_FRF_values_manual,3),'stable');

        % A piston theory supersônica exige U > a_inf (Mach > 1).
        if any(U_FRF_values <= a_inf)
            error(['Todas as velocidades das FRFs devem satisfazer U > a_inf. ' ...
                   'Valor atual de a_inf = %.3f m/s.'],a_inf);
        end

        nU_FRF = numel(U_FRF_values);

        R_FRF_values = R_values;
        nR_FRF = numel(R_FRF_values);
        cmapFRF = lines(nR_FRF);

        legR_FRF = cell(1,nR_FRF);
        for ir = 1:nR_FRF
            legR_FRF{ir} = sprintf('R_L = 10^{%d} \\Omega', ...
                round(log10(R_FRF_values(ir))));
        end

        % ---------------------------------------------------------
        % 9.10.2) Vetor de entrada modal e saida na ponta
        % ---------------------------------------------------------
        Fbase_r = full(PhiA'*Fbase_f);

        tipNodes_FRF = find(abs(nodes(:,1)-L) < 1e-12);
        [~,idxMid_FRF] = min(abs(nodes(tipNodes_FRF,2)-B/2));
        tipNode_FRF = tipNodes_FRF(idxMid_FRF);
        tipDofGlobal_FRF = 3*tipNode_FRF-2;
        [~,tipDofFree_FRF] = ismember(tipDofGlobal_FRF,freeDofs);

        if tipDofFree_FRF == 0
            error('GDL transversal da ponta nao encontrado em freeDofs.');
        end

        phiTip_r = full(PhiA(tipDofFree_FRF,:)); % 1 x nModesAero

        fprintf('\nFRFs piezoaeroelasticas:\n');
        fprintf('Ponta: x = %.6e m | y = %.6e m\n', ...
            nodes(tipNode_FRF,1),nodes(tipNode_FRF,2));
        fprintf('Faixa de frequencia: %.3f a %.3f Hz\n', ...
            freqVec_AERO(1),freqVec_AERO(end));
        fprintf('Velocidades detalhadas: ');
        fprintf('%.3f ',U_FRF_values);
        fprintf('m/s\n');

        % ---------------------------------------------------------
        % 9.10.3) FRFs detalhadas nas velocidades definidas pelo usuario
        % ---------------------------------------------------------
        HrelTip_AERO = complex(NaN(nR_FRF,nFreq_AERO,nU_FRF));
        V_AERO       = complex(NaN(nR_FRF,nFreq_AERO,nU_FRF));
        I_AERO       = complex(NaN(nR_FRF,nFreq_AERO,nU_FRF));
        Pavg_AERO    = NaN(nR_FRF,nFreq_AERO,nU_FRF);

        for iuSel = 1:nU_FRF
            Usel = U_FRF_values(iuSel);

            % Coeficientes da piston theory calculados diretamente para
            % a velocidade escolhida, sem aproximacao pela malha Uvec.
            MachSel = Usel/a_inf;
            if MachSel < 2.0
                warning(['U = %.1f m/s corresponde a Mach %.3f. ' ...
                         'A piston theory empregada deve ser interpretada com cautela abaixo de Mach 2.'], ...
                         Usel,MachSel);
            end
            betaSel = sqrt(MachSel^2-1);
            qinfSel = 0.5*rho_inf*Usel^2;
            lambdaSel = 2*qinfSel/betaSel;
            muSel = (lambdaSel/Usel)* ...
                ((MachSel^2-2)/(MachSel^2-1));

            Ka = lambdaSel*Gxr;
            Ca = muSel*Gvr;
            Kt = Kr + Ka;
            Ct = Cr + Ca;

            fprintf('Calculando FRFs em U = %.3f m/s (Mach %.3f)...\n', ...
                Usel,MachSel);

            modalSummary_FRF = cell(nR_FRF,1);
            stable_FRF = false(nR_FRF,1);

            for ir = 1:nR_FRF
                Rl = R_FRF_values(ir);

                % Verificacao direta da estabilidade na velocidade
                % selecionada. Isso permite usar velocidades fora de Uvec.
                Aaee_FRF = [zeros(nModesAero), eye(nModesAero), zeros(nModesAero,1); ...
                           -(Mr\Kt),          -(Mr\Ct),          Mr\Hr; ...
                            zeros(1,nModesAero), -(Hr'/Cp_total), -1/(Rl*Cp_total)];
                [Vec_FRF,D_FRF] = eig(Aaee_FRF);
                lambda_FRF = diag(D_FRF);
                [sigma_FRF,~,~] = dominant_oscillatory_mode(lambda_FRF);
                modalSummary_FRF{ir} = extract_state_modes(lambda_FRF,Vec_FRF,nModesAero,nModesAero);
                stable_FRF(ir) = sigma_FRF < 0;

                if mascarar_regiao_instavel_FRF && sigma_FRF >= 0
                    fprintf(['  R_L = %.3e Ohm: configuracao instavel; ' ...
                             'FRF estacionaria mascarada.\n'],Rl);
                    continue;
                end

                for iw = 1:nFreq_AERO
                    w = Omega_AERO(iw);
                    Yelec = 1/Rl + 1i*w*Cp_total;

                    Zpiezoaero = -w^2*Mr + 1i*w*Ct + Kt + ...
                        (1i*w/Yelec)*(Hr*Hr');

                    eta = Zpiezoaero\Fbase_r; % resposta por a_b = 1 m/s^2
                    Vp  = -(1i*w/Yelec)*(Hr'*eta);
                    Ip  = Vp/Rl;

                    wTip_acc = phiTip_r*eta;  % w_rel/a_b
                    HrelTip_AERO(ir,iw,iuSel) = w^2*wTip_acc; % w_rel/Y0
                    V_AERO(ir,iw,iuSel) = Vp;
                    I_AERO(ir,iw,iuSel) = Ip;
                    Pavg_AERO(ir,iw,iuSel) = 0.5*abs(Vp)^2/Rl;
                end
            end

            % Figura 2x2 para cada velocidade selecionada.
            figure('Name',sprintf('FRFs piezoaeroelasticas U=%.1f m-s',Usel));
            tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
            title(tl, ...
                sprintf('FRFs piezoaeroelasticas - U_{\\infty} = %.1f m/s',Usel), ...
                'Interpreter','tex');

            nexttile; hold on;
            for ir = 1:nR_FRF
                curve = squeeze(HrelTip_AERO(ir,:,iuSel));
                validCurve = isfinite(curve) & abs(curve) > 0;
                if any(validCurve)
                    semilogx(freqVec_AERO(validCurve), ...
                        20*log10(abs(curve(validCurve))), ...
                        'Color',cmapFRF(ir,:),'LineWidth',1.15);
                end
            end
            grid on; box on;
            xlim([10 2500]);
            xlabel('Frequencia [Hz]');
            ylabel('20log_{10}|w_{rel}(L)/Y_0| [dB]');
            title('Deslocamento relativo da ponta');

            nexttile; hold on;
            for ir = 1:nR_FRF
                curve = squeeze(V_AERO(ir,:,iuSel));
                validCurve = isfinite(curve) & abs(curve) > 0;
                if any(validCurve)
                    semilogx(freqVec_AERO(validCurve), ...
                        20*log10(abs(curve(validCurve))), ...
                        'Color',cmapFRF(ir,:),'LineWidth',1.15);
                end
            end
            grid on; box on;
            xlim([10 2500]);
            xlabel('Frequencia [Hz]');
            ylabel('20log_{10}|V/a_b| [dB]');
            title('Tensao eletrica');

            nexttile; hold on;
            for ir = 1:nR_FRF
                curve = squeeze(I_AERO(ir,:,iuSel));
                validCurve = isfinite(curve) & abs(curve) > 0;
                if any(validCurve)
                    semilogx(freqVec_AERO(validCurve), ...
                        20*log10(abs(curve(validCurve))), ...
                        'Color',cmapFRF(ir,:),'LineWidth',1.15);
                end
            end
            grid on; box on;
            xlim([10 2500]);
            xlabel('Frequencia [Hz]');
            ylabel('20log_{10}|I/a_b| [dB]');
            title('Corrente eletrica');

            nexttile; hold on;
            for ir = 1:nR_FRF
                curve = squeeze(Pavg_AERO(ir,:,iuSel));
                validCurve = isfinite(curve) & curve > 0;
                if any(validCurve)
                    semilogx(freqVec_AERO(validCurve), ...
                        10*log10(curve(validCurve)), ...
                        'Color',cmapFRF(ir,:),'LineWidth',1.15);
                end
            end
            grid on; box on;
            xlim([10 2500]);
            xlabel('Frequencia [Hz]');
            ylabel('10log_{10}(P_{med}/a_b^2) [dB]');
            title('Potencia eletrica media');

            lgd = legend(legR_FRF,'Location','eastoutside','Interpreter','tex');
            lgd.Layout.Tile = 'east';

            if imprimir_modos_console
                fprintf('\n============================================================\n');
                fprintf('MODOS E PICOS DAS FRFs | U = %.3f m/s | Mach = %.3f\n',Usel,MachSel);
                fprintf('============================================================\n');
                for irPrint = 1:nR_FRF
                    Rprint = R_FRF_values(irPrint);
                    modesPrint = modalSummary_FRF{irPrint};
                    if isempty(modesPrint)
                        fprintf('R_L = %.3e Ohm | sem modos oscilatorios validos.\n',Rprint);
                        continue;
                    end
                    [domMode,domFreq,domG,domZeta,domSigma,domPart] = dominant_from_mode_struct(modesPrint);
                    statusTxt = 'ESTAVEL';
                    if ~stable_FRF(irPrint), statusTxt = 'INSTAVEL'; end
                    fprintf(['R_L = %.3e Ohm | %s | modo estrutural dominante = %d ' ...
                             '| f = %.4f Hz | g = %+.4e | zeta = %+.4e ' ...
                             '| Re(lambda) = %+.4e | participacao = %.2f %%\n'], ...
                             Rprint,statusTxt,domMode,domFreq,domG,domZeta,domSigma,100*domPart);

                    if stable_FRF(irPrint)
                        hCurve = squeeze(HrelTip_AERO(irPrint,:,iuSel));
                        vCurve = squeeze(V_AERO(irPrint,:,iuSel));
                        iCurve = squeeze(I_AERO(irPrint,:,iuSel));
                        pCurve = squeeze(Pavg_AERO(irPrint,:,iuSel));
                        print_frf_peaks(freqVec_AERO,hCurve,vCurve,iCurve,pCurve,modesPrint);
                    else
                        fprintf('  FRF estacionaria nao plotada porque a configuracao e instavel.\n');
                    end
                end

                [~,irRef] = min(abs(log10(R_FRF_values)-log10(R_referencia_modos_FRF)));
                if ~isempty(modalSummary_FRF{irRef})
                    labelRef = sprintf('Tabela modal detalhada | U = %.1f m/s | R_L = %.3e Ohm', ...
                        Usel,R_FRF_values(irRef));
                    print_mode_struct(labelRef,modalSummary_FRF{irRef},nModesConsole);
                end
            end
        end

        % ---------------------------------------------------------
        % 9.10.4) Envelopes maximos das FRFs em 1000-2000 m/s
        % ---------------------------------------------------------
        % Esta etapa resume cada FRF pelo seu maior pico na faixa de
        % frequencia. Para limitar o custo computacional, utiliza-se uma
        % malha reduzida de velocidades e frequencias.
        if executar_envelopes_FRF_velocidade

            idxU_env = unique([1:5:nU,nU]);
            Uenv = Uvec(idxU_env);
            nUenv = numel(Uenv);

            freqVec_ENV = logspace(log10(10),log10(2500),550);
            Omega_ENV = 2*pi*freqVec_ENV;

            maxHrel_dB = NaN(nR_FRF,nUenv);
            maxV_dB    = NaN(nR_FRF,nUenv);
            maxI_dB    = NaN(nR_FRF,nUenv);
            maxP_dB    = NaN(nR_FRF,nUenv);

            fprintf('Calculando envelopes maximos das FRFs em 1000-2000 m/s...\n');

            for ku = 1:nUenv
                iuGrid = idxU_env(ku);

                Ka = lambdaA_vec(iuGrid)*Gxr;
                Ca = muA_vec(iuGrid)*Gvr;
                Kt = Kr + Ka;
                Ct = Cr + Ca;

                for ir = 1:nR_FRF
                    Rl = R_FRF_values(ir);

                    if mascarar_regiao_instavel_FRF && ...
                            sigmaMax_AEE(ir,iuGrid) >= 0
                        continue;
                    end

                    peakH = 0;
                    peakV = 0;
                    peakI = 0;
                    peakP = 0;

                    for iw = 1:numel(Omega_ENV)
                        w = Omega_ENV(iw);
                        Yelec = 1/Rl + 1i*w*Cp_total;

                        Zpiezoaero = -w^2*Mr + 1i*w*Ct + Kt + ...
                            (1i*w/Yelec)*(Hr*Hr');

                        eta = Zpiezoaero\Fbase_r;
                        Vp  = -(1i*w/Yelec)*(Hr'*eta);
                        Ip  = Vp/Rl;
                        Hrel = w^2*(phiTip_r*eta);
                        Pmed = 0.5*abs(Vp)^2/Rl;

                        peakH = max(peakH,abs(Hrel));
                        peakV = max(peakV,abs(Vp));
                        peakI = max(peakI,abs(Ip));
                        peakP = max(peakP,Pmed);
                    end

                    maxHrel_dB(ir,ku) = 20*log10(max(peakH,realmin));
                    maxV_dB(ir,ku)    = 20*log10(max(peakV,realmin));
                    maxI_dB(ir,ku)    = 20*log10(max(peakI,realmin));
                    maxP_dB(ir,ku)    = 10*log10(max(peakP,realmin));
                end
            end

            figure('Name','Envelopes das FRFs versus velocidade');
            tl2 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
            title(tl2,'Picos das FRFs piezoaeroelasticas em 1000-2000 m/s');

            nexttile; hold on;
            for ir = 1:nR_FRF
                plot(Uenv,maxHrel_dB(ir,:), ...
                    'Color',cmapFRF(ir,:),'LineWidth',1.25);
            end
            grid on; box on; xlim([1000 2000]);
            xlabel('U_\infty [m/s]');
            ylabel('max 20log_{10}|w_{rel}/Y_0| [dB]');
            title('Pico de deslocamento relativo');

            nexttile; hold on;
            for ir = 1:nR_FRF
                plot(Uenv,maxV_dB(ir,:), ...
                    'Color',cmapFRF(ir,:),'LineWidth',1.25);
            end
            grid on; box on; xlim([1000 2000]);
            xlabel('U_\infty [m/s]');
            ylabel('max 20log_{10}|V/a_b| [dB]');
            title('Pico de tensao');

            nexttile; hold on;
            for ir = 1:nR_FRF
                plot(Uenv,maxI_dB(ir,:), ...
                    'Color',cmapFRF(ir,:),'LineWidth',1.25);
            end
            grid on; box on; xlim([1000 2000]);
            xlabel('U_\infty [m/s]');
            ylabel('max 20log_{10}|I/a_b| [dB]');
            title('Pico de corrente');

            nexttile; hold on;
            for ir = 1:nR_FRF
                plot(Uenv,maxP_dB(ir,:), ...
                    'Color',cmapFRF(ir,:),'LineWidth',1.25);
            end
            grid on; box on; xlim([1000 2000]);
            xlabel('U_\infty [m/s]');
            ylabel('max 10log_{10}(P_{med}/a_b^2) [dB]');
            title('Pico de potencia media');

            lgd2 = legend(legR_FRF,'Location','eastoutside');
            lgd2.Layout.Tile = 'east';
        end
    end

    %% ------------------------------------------------------------
    % 9.11) Observacao para a proxima etapa de harvesting
    % -------------------------------------------------------------
    fprintf('\nIMPORTANTE PARA O HARVESTING:\n');
    fprintf(['A analise linear homogenea determina Ucrit, fcrit e amortecimento, ' ...
             'mas nao determina amplitude absoluta de oscilacao.\n']);
    fprintf(['Para obter tensao e potencia fisicas acima do flutter, a proxima ' ...
             'etapa deve incluir LCO por nao linearidade estrutural/aerodinamica ' ...
             'ou uma excitacao externa especificada.\n']);

end

if executar_FRFs_base_originais

%% ================================================================
% 9) FRFs eletromecanicas para resistencia variavel
% ================================================================

freqVec = linspace(1,1000,5000);
Omega = 2*pi*freqVec;
Rl_values = [1e2, 1e3, 1e4, 1e5, 1e6];
nR = length(Rl_values);
nFreq = length(freqVec);
nFree = length(freeDofs);

tipNodes = find(abs(nodes(:,1) - L) < 1e-12);
[~, idxMid] = min(abs(nodes(tipNodes,2) - B/2));
tipNode = tipNodes(idxMid);
outputDof_global = 3*tipNode - 2;
[~, outputDof_free] = ismember(outputDof_global, freeDofs);
if outputDof_free == 0
    error('O GDL de saida esta preso ou nao existe em freeDofs.');
end

fprintf('\nConfiguracao da saida mecanica:\n');
fprintf('tipNode = %d\n', tipNode);
fprintf('x = %.6e m\n', nodes(tipNode,1));
fprintf('y = %.6e m\n', nodes(tipNode,2));
fprintf('outputDof_global = %d\n', outputDof_global);
fprintf('outputDof_free   = %d\n', outputDof_free);

Ured_elec_FRF  = complex(zeros(nFree, nFreq, nR));
Ufull_elec_FRF = complex(zeros(ndof,  nFreq, nR));

H_disp_acc = complex(zeros(nR,nFreq));
H_trans    = complex(zeros(nR,nFreq));
V_FRF      = complex(zeros(nR,nFreq));
I_FRF      = complex(zeros(nR,nFreq));
Ppeak_FRF  = zeros(nR,nFreq);
Pavg_FRF   = zeros(nR,nFreq);

for r = 1:nR
    Rl = Rl_values(r);
    fprintf('Calculando FRFs para Rl = %.3e Ohm...\n', Rl);

    for i = 1:nFreq
        w = Omega(i);
        Yelec = 1/Rl + 1i*w*Cp_total;

        Zdyn_coupled = -w^2*Mff + 1i*w*Cff + Kff + ...
            (1i*w/Yelec)*(Hff*Hff');

        U = Zdyn_coupled \ Fbase_f;
        Vp = -(1i*w/Yelec)*(Hff'*U);
        Ip = Vp/Rl;

        Ured_elec_FRF(:,i,r) = U;
        Ufull_elec_FRF(freeDofs,i,r) = U;

        H_disp_acc(r,i) = U(outputDof_free);
        H_trans(r,i) = (w^2)*U(outputDof_free);
        V_FRF(r,i) = Vp;
        I_FRF(r,i) = Ip;
        Ppeak_FRF(r,i) = abs(Vp)^2/Rl;
        Pavg_FRF(r,i)  = 0.5*Ppeak_FRF(r,i);
    end
end

Wnodes_elec_FRF = Ufull_elec_FRF(1:3:end,:,:);

%% ================================================================
% 10) FRF mecanica desacoplada para comparacao
% ================================================================

Ured_unc_FRF  = complex(zeros(nFree,nFreq));
Ufull_unc_FRF = complex(zeros(ndof,nFreq));
H_disp_acc_unc = complex(zeros(1,nFreq));
H_trans_unc    = complex(zeros(1,nFreq));

for i = 1:nFreq
    w = Omega(i);
    Zdyn_unc = -w^2*Mff + 1i*w*Cff + Kff;
    Uunc = Zdyn_unc \ Fbase_f;

    Ured_unc_FRF(:,i) = Uunc;
    Ufull_unc_FRF(freeDofs,i) = Uunc;
    H_disp_acc_unc(i) = Uunc(outputDof_free);
    H_trans_unc(i) = w^2*Uunc(outputDof_free);
end

Wnodes_unc_FRF = Ufull_unc_FRF(1:3:end,:);

%% ================================================================
% 11) Plots das FRFs em dB
% ================================================================
% Conversoes usadas:
%
% Grandezas de amplitude, como tensao, corrente e deslocamento:
%   X_dB = 20*log10(|X|)
%
% Grandezas de potencia:
%   P_dB = 10*log10(P)
%
% Referencias implicitas:
% - tensao:            dB re 1 [V/(m/s^2)]
% - corrente:          dB re 1 [A/(m/s^2)]
% - potencia:          dB re 1 [W/(m/s^2)^2]
% - deslocamento:      dB re 1 [m/(m/s^2)]
% - transmissibilidade dB re 1 [-]

legendas_R = cell(1,nR);
for r = 1:nR
    legendas_R{r} = sprintf('R_L = 10^{%d} \\Omega', round(log10(Rl_values(r))));
end

legendas_R_unc = [legendas_R, {'Uncoupled'}];

% ----------------------------------------------------------------
% 11.1) Tensao eletrica em dB
% ----------------------------------------------------------------

figure; hold on;
for r = 1:nR
    V_dB = 20*log10(max(abs(V_FRF(r,:)), realmin));
    plot(freqVec, V_dB, 'LineWidth', 1.2);
end
grid on;
xlabel('Frequencia [Hz]');
ylabel('Magnitude [dB re 1 V/(m/s^2)]');
title('FRF de tensao eletrica em dB');
legend(legendas_R, 'Location', 'best');
xlim([0 1000]);

% ----------------------------------------------------------------
% 11.2) Corrente eletrica em dB
% ----------------------------------------------------------------

figure; hold on;
for r = 1:nR
    I_dB = 20*log10(max(abs(I_FRF(r,:)), realmin));
    plot(freqVec, I_dB, 'LineWidth', 1.2);
end
grid on;
xlabel('Frequencia [Hz]');
ylabel('Magnitude [dB re 1 A/(m/s^2)]');
title('FRF de corrente eletrica em dB');
legend(legendas_R, 'Location', 'best');
xlim([0 1000]);

% ----------------------------------------------------------------
% 11.3) Potencia eletrica de pico em dB
% ----------------------------------------------------------------

figure; hold on;
for r = 1:nR
    Ppeak_dB = 10*log10(max(Ppeak_FRF(r,:), realmin));
    plot(freqVec, Ppeak_dB, 'LineWidth', 1.2);
end
grid on;
xlabel('Frequencia [Hz]');
ylabel('Magnitude [dB re 1 W/(m/s^2)^2]');
title('FRF de potencia eletrica de pico em dB');
legend(legendas_R, 'Location', 'best');
xlim([0 1000]);

% ----------------------------------------------------------------
% 11.4) Potencia eletrica media em dB
% ----------------------------------------------------------------

figure; hold on;
for r = 1:nR
    Pavg_dB = 10*log10(max(Pavg_FRF(r,:), realmin));
    plot(freqVec, Pavg_dB, 'LineWidth', 1.2);
end
grid on;
xlabel('Frequencia [Hz]');
ylabel('Magnitude [dB re 1 W/(m/s^2)^2]');
title('FRF de potencia eletrica media em dB');
legend(legendas_R, 'Location', 'best');
xlim([0 1000]);

% ----------------------------------------------------------------
% 11.5) Deslocamento da ponta por aceleracao de base em dB
% ----------------------------------------------------------------

figure; hold on;
for r = 1:nR
    Hdisp_dB = 20*log10(max(abs(H_disp_acc(r,:)), realmin));
    plot(freqVec, Hdisp_dB, 'LineWidth', 1.2);
end
Hdisp_unc_dB = 20*log10(max(abs(H_disp_acc_unc), realmin));
plot(freqVec, Hdisp_unc_dB, 'k--', 'LineWidth', 1.2);
grid on;
xlabel('Frequencia [Hz]');
ylabel('Magnitude [dB re 1 m/(m/s^2)]');
title('FRF de deslocamento por aceleracao de base em dB');
legend(legendas_R_unc, 'Location', 'best');
xlim([0 1000]);

% ----------------------------------------------------------------
% 11.6) Transmissibilidade mecanica da ponta em dB
% ----------------------------------------------------------------

figure; hold on;
for r = 1:nR
    Htrans_dB = 20*log10(max(abs(H_trans(r,:)), realmin));
    plot(freqVec, Htrans_dB, 'LineWidth', 1.2);
end
Htrans_unc_dB = 20*log10(max(abs(H_trans_unc), realmin));
plot(freqVec, Htrans_unc_dB, 'k--', 'LineWidth', 1.2);
grid on;
xlabel('Frequencia [Hz]');
ylabel('Magnitude [dB re 1]');
title('Transmissibilidade mecanica da ponta em dB');
legend(legendas_R_unc, 'Location', 'best');
xlim([0 1000]);

% ----------------------------------------------------------------
% 11.7) Zoom da transmissibilidade no primeiro modo em dB
% ----------------------------------------------------------------

figure; hold on;
for r = 1:nR
    Htrans_dB = 20*log10(max(abs(H_trans(r,:)), realmin));
    plot(freqVec, Htrans_dB, 'LineWidth', 1.2);
end
plot(freqVec, Htrans_unc_dB, 'k--', 'LineWidth', 1.2);
grid on;
xlabel('Frequencia [Hz]');
ylabel('Magnitude [dB re 1]');
title('Zoom da transmissibilidade no primeiro modo em dB');
legend(legendas_R_unc, 'Location', 'best');
xlim([45 52]);

% ----------------------------------------------------------------
% 11.8) Zoom da potencia eletrica de pico no primeiro modo em dB
% ----------------------------------------------------------------

figure; hold on;
for r = 1:nR
    Ppeak_dB = 10*log10(max(Ppeak_FRF(r,:), realmin));
    plot(freqVec, Ppeak_dB, 'LineWidth', 1.2);
end
grid on;
xlabel('Frequencia [Hz]');
ylabel('Magnitude [dB re 1 W/(m/s^2)^2]');
title('Zoom da potencia eletrica de pico no primeiro modo em dB');
legend(legendas_R, 'Location', 'best');
xlim([45 52]);

%% ================================================================
% 12) Campo de deslocamento em dB para uma frequencia e resistencia escolhidas
% ================================================================

[~, idxModo1] = min(abs(freqVec - freq_sc(1)));
[~, idxR_plot] = min(abs(Rl_values - 1e4));

w_field = Wnodes_elec_FRF(:, idxModo1, idxR_plot);
Wplot_FRF_dB = reshape(20*log10(max(abs(w_field), realmin)), [nnx,nny])';

figure;
surf(xGrid, yGrid, Wplot_FRF_dB);
shading interp;
xlabel('x [m]');
ylabel('y [m]');
zlabel('Magnitude [dB re 1 m/(m/s^2)]');
title(sprintf('Campo de deslocamento em dB: f = %.2f Hz, R_L = %.0e Ohm', ...
    freqVec(idxModo1), Rl_values(idxR_plot)));
view(35,25);
grid on;
colorbar;

end % executar_FRFs_base_originais

%% =================================================================
% FUNCOES LOCAIS
% =================================================================

function [Gxe, Gve] = kirchhoff_piston_element(a, b, ngp)
    % Matrizes geometricas elementares para piston theory com escoamento +x.
    % Gxe e geralmente nao simetrica e NAO deve ser simetrizada.
    Gxe = zeros(12,12);
    Gve = zeros(12,12);

    [gp, gw] = gauss_points(ngp);
    A = acm_A_matrix(a,b);
    J = a*b/4;

    for i = 1:ngp
        for j = 1:ngp
            xi = gp(i);
            eta = gp(j);
            weight = gw(i)*gw(j);

            x = (a/2)*(xi+1);
            y = (b/2)*(eta+1);

            [Nw, Nwx, ~, ~, ~, ~] = acm_shape_functions(x,y,A);

            Gxe = Gxe + weight*J*(Nw'*Nwx);
            Gve = Gve + weight*J*(Nw'*Nw);
        end
    end
end

function [sigmaMax, freqDom, gDom] = dominant_oscillatory_mode(lambda)
    % Seleciona, entre as raizes com parte imaginaria positiva, aquela com
    % maior parte real. Esta raiz define o envelope de estabilidade.
    tol = 1e-8;
    idx = find(imag(lambda) > tol & isfinite(real(lambda)) & isfinite(imag(lambda)));

    if isempty(idx)
        sigmaMax = NaN;
        freqDom = NaN;
        gDom = NaN;
        return;
    end

    lamOsc = lambda(idx);
    [sigmaMax,iloc] = max(real(lamOsc));
    lamDom = lamOsc(iloc);
    freqDom = abs(imag(lamDom))/(2*pi);
    gDom = 2*real(lamDom)/max(abs(imag(lamDom)),eps);
end

function [lamSel, vecSel] = track_oscillatory_modes(lambda, Vec, nMech, nTrack, prevLam, prevVec)
    % Rastreamento modal por MAC dos deslocamentos mecanicos e proximidade
    % do autovalor. Evita simples trocas de ordem perto da coalescencia.
    tol = 1e-8;
    idx = find(imag(lambda) > tol & isfinite(real(lambda)) & isfinite(imag(lambda)));

    if isempty(idx)
        lamSel = [];
        vecSel = [];
        return;
    end

    lamCand = lambda(idx);
    vecCand = Vec(:,idx);

    [~,ord] = sort(abs(imag(lamCand)),'ascend');
    lamCand = lamCand(ord);
    vecCand = vecCand(:,ord);

    nCand = min(numel(lamCand), max(4*nTrack,nTrack));
    lamCand = lamCand(1:nCand);
    vecCand = vecCand(:,1:nCand);

    % Normaliza cada autovetor pela parte mecanica de deslocamentos.
    for j = 1:nCand
        nq = norm(vecCand(1:nMech,j));
        if nq > eps
            vecCand(:,j) = vecCand(:,j)/nq;
        end
    end

    nOut = min(nTrack,nCand);

    if isempty(prevLam) || isempty(prevVec)
        lamSel = lamCand(1:nOut);
        vecSel = vecCand(:,1:nOut);
        return;
    end

    nPrev = min([numel(prevLam),size(prevVec,2),nOut]);
    lamSel = NaN(nPrev,1);
    vecSel = complex(NaN(size(Vec,1),nPrev));
    used = false(1,nCand);

    for i = 1:nPrev
        qPrev = prevVec(1:nMech,i);
        scores = -inf(1,nCand);

        for j = 1:nCand
            if used(j)
                continue;
            end

            qCur = vecCand(1:nMech,j);
            den = real((qPrev'*qPrev)*(qCur'*qCur));
            if den <= eps
                mac = 0;
            else
                mac = abs(qPrev'*qCur)^2/den;
            end

            freqScale = max(abs(imag(prevLam(i))),1);
            df = abs(imag(lamCand(j))-imag(prevLam(i)))/freqScale;
            ds = abs(real(lamCand(j))-real(prevLam(i)))/freqScale;

            scores(j) = mac/(1 + 0.25*df + 0.05*ds);
        end

        [~,jBest] = max(scores);
        used(jBest) = true;
        lamSel(i) = lamCand(jBest);
        vecSel(:,i) = vecCand(:,jBest);
    end
end

function [Ucrit, MachCrit, fcrit, idxCrit] = find_flutter_crossing(Uvec, sigmaMax, freqDom, aInf)
    % Interpolacao linear do primeiro cruzamento de Re(lambda) de negativo
    % para nao negativo.
    Ucrit = NaN;
    MachCrit = NaN;
    fcrit = NaN;
    idxCrit = NaN;

    valid = isfinite(sigmaMax) & isfinite(Uvec);
    if nnz(valid) < 2
        return;
    end

    if sigmaMax(1) >= 0
        Ucrit = Uvec(1);
        MachCrit = Ucrit/aInf;
        fcrit = freqDom(1);
        idxCrit = 1;
        return;
    end

    idx = find(sigmaMax(1:end-1) < 0 & sigmaMax(2:end) >= 0,1,'first');
    if isempty(idx)
        return;
    end

    s1 = sigmaMax(idx);
    s2 = sigmaMax(idx+1);
    U1 = Uvec(idx);
    U2 = Uvec(idx+1);

    if abs(s2-s1) <= eps
        alpha = 0;
    else
        alpha = -s1/(s2-s1);
    end
    alpha = min(max(alpha,0),1);

    Ucrit = U1 + alpha*(U2-U1);
    MachCrit = Ucrit/aInf;

    if isfinite(freqDom(idx)) && isfinite(freqDom(idx+1))
        fcrit = freqDom(idx) + alpha*(freqDom(idx+1)-freqDom(idx));
    end
    idxCrit = idx;
end

function [Ke, Me] = kirchhoff_rect_element(a, b, layers, ngp)
    Ke = zeros(12,12);
    Me = zeros(12,12);
    [gp, gw] = gauss_points(ngp);
    A = acm_A_matrix(a,b);
    J = a*b/4;
    for i = 1:ngp
        for j = 1:ngp
            xi  = gp(i);
            eta = gp(j);
            weight = gw(i)*gw(j);
            x = (a/2)*(xi + 1);
            y = (b/2)*(eta + 1);
            [C, Cx, Cy, Cxx, Cyy, Cxy] = acm_shape_functions(x,y,A);
            Bk = [Cxx; Cyy; 2*Cxy];
            for k = 1:length(layers)
                z_bot = layers(k).z_bot;
                z_top = layers(k).z_top;
                Q = layers(k).Q;
                rho = layers(k).rho;
                I0 = z_top - z_bot;
                I2 = (z_top^3 - z_bot^3)/3;
                Db = I2*Q;
                Ke_layer = Bk'*Db*Bk;
                Me_layer = rho*(I0*(C'*C) + I2*(Cx'*Cx + Cy'*Cy));
                Ke = Ke + weight*J*Ke_layer;
                Me = Me + weight*J*Me_layer;
            end
        end
    end
end

function Fe = kirchhoff_base_force_element(a, b, layers, ngp)
    Fe = zeros(12,1);
    [gp, gw] = gauss_points(ngp);
    A = acm_A_matrix(a,b);
    J = a*b/4;
    m_area = 0;
    for k = 1:length(layers)
        h_layer = layers(k).z_top - layers(k).z_bot;
        m_area = m_area + layers(k).rho*h_layer;
    end
    for i = 1:ngp
        for j = 1:ngp
            xi = gp(i);
            eta = gp(j);
            weight = gw(i)*gw(j);
            x = (a/2)*(xi + 1);
            y = (b/2)*(eta + 1);
            [C, ~, ~, ~, ~, ~] = acm_shape_functions(x,y,A);
            Fe = Fe + weight*J*C'*m_area;
        end
    end
end

function [He, Cpe] = kirchhoff_piezo_element(a, b, layers, ngp)
    He = zeros(12,1);
    Cpe = 0;
    [gp, gw] = gauss_points(ngp);
    A = acm_A_matrix(a,b);
    J = a*b/4;
    for i = 1:ngp
        for j = 1:ngp
            xi = gp(i);
            eta = gp(j);
            weight = gw(i)*gw(j);
            x = (a/2)*(xi + 1);
            y = (b/2)*(eta + 1);
            [~, ~, ~, Cxx, Cyy, ~] = acm_shape_functions(x,y,A);
            for k = 1:length(layers)
                if ~layers(k).isPiezo
                    continue;
                end
                z_bot = layers(k).z_bot;
                z_top = layers(k).z_top;
                hp_layer = z_top - z_bot;
                e31bar = layers(k).e31bar;
                e32bar = layers(k).e32bar;
                eps33bar = layers(k).eps33bar;
                I1 = (z_top^2 - z_bot^2)/2;
                He_local = (I1/hp_layer)*(e31bar*Cxx' + e32bar*Cyy');
                Cpe_local = eps33bar/hp_layer;
                He = He + weight*J*He_local;
                Cpe = Cpe + weight*J*Cpe_local;
            end
        end
    end
end

function A = acm_A_matrix(a,b)
    nodeCoords = [0 0; a 0; a b; 0 b];
    A = zeros(12,12);
    row = 0;
    for k = 1:4
        x = nodeCoords(k,1);
        y = nodeCoords(k,2);
        P  = poly_P(x,y);
        Px = poly_dPdx(x,y);
        Py = poly_dPdy(x,y);
        row = row + 1; A(row,:) = P;
        row = row + 1; A(row,:) = Py;
        row = row + 1; A(row,:) = Px;
    end
end

function [C, Cx, Cy, Cxx, Cyy, Cxy] = acm_shape_functions(x,y,A)
    P   = poly_P(x,y);
    Px  = poly_dPdx(x,y);
    Py  = poly_dPdy(x,y);
    Pxx = poly_d2Pdx2(x,y);
    Pyy = poly_d2Pdy2(x,y);
    Pxy = poly_d2Pdxdy(x,y);
    C   = P/A;
    Cx  = Px/A;
    Cy  = Py/A;
    Cxx = Pxx/A;
    Cyy = Pyy/A;
    Cxy = Pxy/A;
end

function P = poly_P(x,y)
    P = [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, x^3*y, x*y^3];
end

function Px = poly_dPdx(x,y)
    Px = [0, 1, 0, 2*x, y, 0, 3*x^2, 2*x*y, y^2, 0, 3*x^2*y, y^3];
end

function Py = poly_dPdy(x,y)
    Py = [0, 0, 1, 0, x, 2*y, 0, x^2, 2*x*y, 3*y^2, x^3, 3*x*y^2];
end

function Pxx = poly_d2Pdx2(x,y)
    Pxx = [0, 0, 0, 2, 0, 0, 6*x, 2*y, 0, 0, 6*x*y, 0];
end

function Pyy = poly_d2Pdy2(x,y)
    Pyy = [0, 0, 0, 0, 0, 2, 0, 0, 2*x, 6*y, 0, 6*x*y];
end

function Pxy = poly_d2Pdxdy(x,y)
    Pxy = [0, 0, 0, 0, 1, 0, 0, 2*x, 2*y, 0, 3*x^2, 3*y^2];
end

function [gp, gw] = gauss_points(n)
    switch n
        case 2
            gp = [-1/sqrt(3), 1/sqrt(3)];
            gw = [1, 1];
        case 3
            gp = [-sqrt(3/5), 0, sqrt(3/5)];
            gw = [5/9, 8/9, 5/9];
        case 4
            gp = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];
            gw = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454];
        case 5
            gp = [-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664];
            gw = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189];
        otherwise
            error('Numero de pontos de Gauss nao implementado.');
    end
end



function modes = extract_state_modes(lambda,Vec,nMech,nOut)
    % Extrai os ramos oscilatorios com parte imaginaria positiva e associa
    % cada ramo ao modo estrutural da base PhiA com maior participacao.
    tol = 1e-8;
    idx = find(imag(lambda) > tol & isfinite(real(lambda)) & isfinite(imag(lambda)));

    modes = struct('branch',{},'structuralMode',{},'participation',{}, ...
                   'frequency',{},'sigma',{},'g',{},'zeta',{},'lambda',{});
    if isempty(idx)
        return;
    end

    lam = lambda(idx);
    vec = Vec(:,idx);
    [~,ord] = sort(abs(imag(lam)),'ascend');
    lam = lam(ord);
    vec = vec(:,ord);
    nKeep = min(nOut,numel(lam));

    for k = 1:nKeep
        q = vec(1:nMech,k);
        energy = abs(q).^2;
        den = sum(energy);
        if den <= eps
            structuralMode = NaN;
            participation = NaN;
        else
            [participation,structuralMode] = max(energy/den);
        end
        lk = lam(k);
        modes(k).branch = k;
        modes(k).structuralMode = structuralMode;
        modes(k).participation = participation;
        modes(k).frequency = abs(imag(lk))/(2*pi);
        modes(k).sigma = real(lk);
        modes(k).g = 2*real(lk)/max(abs(imag(lk)),eps);
        modes(k).zeta = -real(lk)/max(abs(lk),eps);
        modes(k).lambda = lk;
    end
end

function print_state_modes(label,Astate,nMech,nOut)
    [Vec,D] = eig(Astate);
    modes = extract_state_modes(diag(D),Vec,nMech,nOut);
    print_mode_struct(label,modes,nOut);
end

function print_mode_struct(label,modes,nOut)
    fprintf('\n%s\n',label);
    fprintf('Ramo  Modo estrutural  Participacao[%%]   f[Hz]       Re(lambda)        g             zeta\n');
    fprintf('----  ---------------  ---------------  ----------  ---------------  ------------  ------------\n');
    nPrint = min(nOut,numel(modes));
    for k = 1:nPrint
        fprintf('%4d  %15d  %15.3f  %10.4f  %+15.6e  %+12.5e  %+12.5e\n', ...
            modes(k).branch,modes(k).structuralMode,100*modes(k).participation, ...
            modes(k).frequency,modes(k).sigma,modes(k).g,modes(k).zeta);
    end
end

function print_dominant_state_mode(label,Astate,nMech,Uvalue)
    [Vec,D] = eig(Astate);
    modes = extract_state_modes(diag(D),Vec,nMech,10*nMech);
    if isempty(modes)
        fprintf('%s | U = %.4f m/s | nenhum modo oscilatorio valido.\n',label,Uvalue);
        return;
    end
    [modeNum,freq,g,zeta,sigma,part] = dominant_from_mode_struct(modes);
    fprintf(['%s | Ucrit = %.4f m/s | modo estrutural dominante = %d ' ...
             '| f = %.4f Hz | Re(lambda) = %+.5e | g = %+.5e ' ...
             '| zeta = %+.5e | participacao = %.2f %%\n'], ...
             label,Uvalue,modeNum,freq,sigma,g,zeta,100*part);
end

function [modeNum,freq,g,zeta,sigma,part] = dominant_from_mode_struct(modes)
    sigmaAll = [modes.sigma];
    [~,idx] = max(sigmaAll);
    modeNum = modes(idx).structuralMode;
    freq = modes(idx).frequency;
    g = modes(idx).g;
    zeta = modes(idx).zeta;
    sigma = modes(idx).sigma;
    part = modes(idx).participation;
end

function print_frf_peaks(freqVec,hCurve,vCurve,iCurve,pCurve,modes)
    [fH,valH] = max_valid_peak(freqVec,abs(hCurve));
    [fV,valV] = max_valid_peak(freqVec,abs(vCurve));
    [fI,valI] = max_valid_peak(freqVec,abs(iCurve));
    [fP,valP] = max_valid_peak(freqVec,pCurve);

    mH = nearest_structural_mode(fH,modes);
    mV = nearest_structural_mode(fV,modes);
    mI = nearest_structural_mode(fI,modes);
    mP = nearest_structural_mode(fP,modes);

    fprintf(['  Picos: desloc. f=%.3f Hz (modo %d, %.3e) | ' ...
             'tensao f=%.3f Hz (modo %d, %.3e) | ' ...
             'corrente f=%.3f Hz (modo %d, %.3e) | ' ...
             'potencia f=%.3f Hz (modo %d, %.3e)\n'], ...
             fH,mH,valH,fV,mV,valV,fI,mI,valI,fP,mP,valP);
end

function [fPeak,valPeak] = max_valid_peak(freqVec,values)
    valid = isfinite(values) & values > 0;
    if ~any(valid)
        fPeak = NaN;
        valPeak = NaN;
        return;
    end
    idxValid = find(valid);
    [valPeak,iLocal] = max(values(valid));
    idxPeak = idxValid(iLocal);
    fPeak = freqVec(idxPeak);
end

function modeNum = nearest_structural_mode(fPeak,modes)
    if ~isfinite(fPeak) || isempty(modes)
        modeNum = NaN;
        return;
    end
    frequencies = [modes.frequency];
    [~,idx] = min(abs(frequencies-fPeak));
    modeNum = modes(idx).structuralMode;
end

function [zNeutral, n, hpa, hsa, hpc, ha, hb, hc] = compute_neutral_axis_transformed_section(hs, hp, Ys, Yp)
    n = Ys/Yp;
    hpa = (hp^2 + 2*n*hp*hs + n*hs^2)/(2*(hp + n*hs));
    hsa = (hp^2 + 2*hp*hs + n*hs^2)/(2*(hp + n*hs));
    hpc = (n*hs*(hp + hs))/(2*(hp + n*hs));
    ha = -hsa;
    hb = hpa - hp;
    hc = hpa;
    zNeutral = hsa;
end
