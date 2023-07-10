function TranscriptionTranslation

clear all
close all hidden
warning off

%Set lognormal functions and noise characteristics
CV = 0.1;

sigma = sqrt(log(CV^2+1));

numsims = 100;
numtps = 10;
rng(10000)
KTXN = 1;
KTRL = 1;
RDEG = 0.1;
PDEG = 0.01;
for i = 1:numsims
    %Set rate parameters
    k_txn = lognrnd(muFunc(KTXN,(CV*KTXN)^2),sigma);
    k_trl = lognrnd(muFunc(KTRL,(CV*KTRL)^2),sigma);
    k_RNAdeg = lognrnd(muFunc(RDEG,(CV*RDEG)^2),sigma);
    k_protdeg = lognrnd(muFunc(PDEG,(CV*PDEG)^2),sigma);
    Model_params = [k_txn k_trl k_RNAdeg k_protdeg];

    %Set initial conditions at 50% the steady-state values
    RNA_IC = 0.5*k_txn/k_RNAdeg;
    Prot_IC = 0.5*RNA_IC*k_trl/k_protdeg;
    Gene_IC = [RNA_IC Prot_IC];

    %Set time interval
    tspan = 100;

    [Time,Gene_soln] = ode15s(@(t,y) GeneODE(t,y,Model_params),...
        linspace(0,tspan,tspan+1), Gene_IC,[]);
    Register(:,:,i)=Gene_soln(randi(tspan+1,1,numtps),:);
end

Register=reshape(shiftdim(Register,2),numsims*numtps,size(Register,2));
x = Register(:,1);
y = Register(:,2);
plot(log10(x),log10(y),'b.')
hold on

hyperlin = @(a,xdata) a(1).*(a(2).*xdata./(xdata + a(3)) + xdata);
a0 = [1,1,1];

% First, fit with lsqcurvefit since that can restrict C > 0. This also
% provides a better initial guess to nlinfit, which can robustly weight the
% cost function when calculating the error for the fit.
options_hl = optimoptions('lsqcurvefit','FiniteDifferenceType','central','Display','off');
beta0_hl = lsqcurvefit(hyperlin, a0, x, y, [-Inf,-Inf,0], [], options_hl);

% nlinfit is the better fitter for calculating error, because a logistic
% weighting fucntion can be set.
opts = statset('nlinfit');
opts.RobustWgtFun = 'logistic';
opts.Display = 'off';
mdl_hl = fitnlm(x, y, hyperlin, beta0_hl, 'Options', opts);

% The problem with nlinfit is that C > 0 cannot be enforced. In the case
% that this does occur, use the output from nlinfit as the initial guess
% for, once again, lsqcurvefit to restrict C > 0.
if (mdl_hl.Coefficients.Estimate(3) < 0)
    beta = lsqcurvefit(hyperlin, mdl_hl.Coefficients.Estimate(:), x, y, [-Inf -Inf 0], [], options_hl);
else
    beta = mdl_hl.Coefficients.Estimate(:);
end
% subplot(1,2,1)
% plot(sort(x), hyperlin(beta, sort(x)), 'g-')
% subplot(1,2,2)
plot(sort(log10(x)), log10(hyperlin(beta,sort(x))), 'g-')
xlim([0.5 1.5])
ylim([2 3])
xlabel('RNA')
ylabel('Protein')

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dGene_dt = GeneODE(t,Gene,params)

param_cell = num2cell(params);
[k_r,k_p,k_rd,k_pd] = param_cell{:};

Gene_cell = num2cell(Gene);
[RNA,Protein] = Gene_cell{:};

dRNA_dt = k_r - k_rd*RNA;
dProtein_dt = k_p*RNA - k_pd*Protein;

dGene_dt = [dRNA_dt; dProtein_dt];

return

function mu = muFunc(M,V)
    mu = log((M^2)/sqrt(V+M^2));
return
