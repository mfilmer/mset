% SET.Cs = 30e-18;
% SET.Cd = 30e-18;
% SET.Cg = 0.1e-18;
% SET.Gd = 1e-6;
% SET.Gs = 1e-6;
% SET.T = 0.3;
% Bias.Vs = 0;
% Bias.Vd = 0;
% Bias.Vg = 0;
function [G, vds, vgs] = basicset(SET, Bias)
    nvg = 201;
    vgs = linspace(-3000e-3,3000e-3,nvg);
    nvd = 201;
    % Do some math so that the start and end points after a numeric
    % derivative become the desired start and end points
    vds_min = -20e-3;
    vds_max = -vds_min;
    vds_step = (vds_max - vds_min)/(nvd-1);
    vds = linspace(vds_min-vds_step/2,vds_max+vds_step/2,nvd+1);
    
    I = zeros(nvd,nvg);
    %tCount = zeros(nvd,nvg);
    
    %Bias.Vs = 0;
    for ivg = 1:nvg
        h = waitbar(ivg/nvg);
        Bias.Vg = vgs(ivg);
        for ivd = 1:nvd+1
            Bias.Vd = vds(ivd);
            
            %Noptions = -5:5;
            %for Ninitial = -5:5
                %Nvalid = findTransitions(SET, Bias, Ninitial, Noptions);
                %checkReturnTransition = @(N) isValidTransition(SET, Bias, N, Ninitial);
                %validCount = sum(arrayfun(checkReturnTransition,Nvalid));
                %tCount(ivd,ivg) = tCount(ivd,ivg) + validCount;
            %end
            I(ivd,ivg) = current(SET, Bias);
        end
    end
    
    close(h);
    
    G = diff(I,1) ./ repmat(diff(vds)',1,nvg);
    vds = linspace(vds_min, vds_max, nvd);
    %G = tCount./2;
end

function Nvalid = findTransitions(SET, Bias, Ninitial, Noptions)
    checkTransition = @(N) isValidTransition(SET, Bias, Ninitial, N);
    Nvalid = Noptions(arrayfun(checkTransition, Noptions));
end

% Returns true if the transition has a negative energy cost
function isValid = isValidTransition(SET, Bias, Ninitial, Nfinal)
    e = 1.60217e-19;            % Coulombs
    % It isn't a transition if nothing changes
    if Ninitial == Nfinal
        isValid = false;
        return;
    end
    % We only deal with single step transitions here
    if abs(Ninitial-Nfinal) > 1
        isValid = false;
        return
        oneStep = sign(Nfinal-Ninitial);
        isValid = isValidTransition(SET, Bias, Ninitial, Ninitial+oneStep) && ...
                isValidTransition(SET, Bias, Ninitial+oneStep, Nfinal);
        return
    end
    dN = Nfinal - Ninitial;
    dE_source = E_el(SET, Bias, Nfinal) - E_el(SET, Bias, Ninitial) + dN*e*Bias.Vs;
    dE_drain = E_el(SET, Bias, Nfinal) - E_el(SET, Bias, Ninitial) + dN*e*Bias.Vd;
    if dE_source < 0 || dE_drain < 0
        isValid = true;
    else
        isValid = false;
    end
end

% Equ 3.13
function E = E_el(SET, Bias, N)
    e = 1.60217e-19;            % Coulombs
    Ctotal = SET.Cs + SET.Cd + SET.Cg;
    Ec = e^2/(2*Ctotal);
    q = SET.Cs * Bias.Vs + SET.Cd * Bias.Vd + SET.Cg * Bias.Vg;
    %E = Ec.*(N - q/e).^2;
    E = Ec.*(N - q/e).^2 + SET.DeltaI * mod(N,2);
end

function gamma = gamma_L(SET, Bias, N, dN)
    e = 1.60217e-19;            % Coulombs
    kT = 1.38065e-23 * SET.T;   % Joules
    
    if abs(dN) ~= 1
        gamma = 0;
        return;
    end
    
    % Based on Equ 3.14
    dE = E_el(SET, Bias, N+dN) - E_el(SET, Bias, N) + dN*e*Bias.Vs;
    % Add Superconducting Delta
    if dN > 1
        dE = dE + 2*SET.DeltaI;
    elseif dN < 1
        dE = dE + 2*SET.DeltaL;
    end
    
    if dE == 0
        gamma = SET.Gs./e.^2;
        return;
    end
    gamma = SET.Gs./e.^2 .* dE ./ (exp(dE/kT) - 1);
end

function gamma = gamma_R(SET, Bias, N, dN)
    e = 1.60217e-19;            % Coulombs
    kT = 1.38065e-23 * SET.T;   % Joules
    
    if abs(dN) ~= 1
        gamma = 0;
        return;
    end
    
    % Based on Equ 3.14
    dE = E_el(SET, Bias, N+dN) - E_el(SET, Bias, N) + dN*e*Bias.Vd;
    % Add Superconducting Delta
    if dN > 1
        dE = dE + 2*SET.DeltaI;
    elseif dN < 1
        dE = dE + 2*SET.DeltaL;
    end
    
    if dE == 0
        gamma = SET.Gd./e.^2;
        return;
    end
    gamma = SET.Gd./e.^2 .* dE ./ (exp(dE/kT) - 1);
end

% Equ 3.20
function I = current(SET, Bias)
    e = 1.60217e-19;            % Coulombs
    
    % These should be equal, and it seems that they are
    %gammaDiff = @(N) sum(gamma_L(SET, Bias, N, 1)) - sum(gamma_L(SET, Bias, N, -1));
    gammaDiff = @(N) sum(gamma_R(SET, Bias, N, -1)) - sum(gamma_R(SET, Bias, N, 1));
    
    gamma = @(N, dN) gamma_L(SET, Bias, N, dN) + gamma_R(SET, Bias, N, dN);
    
    % Adapted from Equ 3.46
    %gammaMat = [0            gamma(-1, -1)  0            0            0;
    %            gamma(-2,1)  0              gamma(0,-1)  0            0;
    %            0            gamma(-1,1)    0            gamma(1,-1)  0;
    %            0            0              gamma(0,1)   0            gamma(2,-1);
    %            0            0              0            gamma(1,1)   0];
    ns = -5:5;
    gammaPlus = @(N) gamma(N, 1);
    gammaMinus = @(N) gamma(N, -1);
    gammaMat = zeros(length(ns));
    gammaMat = gammaMat + diag(arrayfun(gammaPlus, ns(1:end-1)), -1);
    gammaMat = gammaMat + diag(arrayfun(gammaMinus, ns(2:end)), 1);
    gammaMat = gammaMat - diag(sum(gammaMat,1));
    % Equ 3.22
    [V,D] = eig(gammaMat);
    D = diag(D);
    [~,index] = min(abs(D));
    pns = V(:,index);
    pns = pns/sum(pns);
    
    I = e * sum(arrayfun(gammaDiff, ns)'.*pns);
end
