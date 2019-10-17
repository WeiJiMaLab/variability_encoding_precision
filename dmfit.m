function fitpars = dmfit(x)

% find initial parameter estimate
maxLLH=-Inf;
wG_vec = linspace(0,1,10);
kappa_vec = linspace(1,20,10);

for ii = randperm(length(wG_vec))
    for jj = randperm(length(kappa_vec))
        wG_hat = wG_vec(ii);
        kappa_hat = kappa_vec(jj);
        LLH = -mLLHfun([wG_hat kappa_hat],x);
        if LLH>maxLLH
            maxLLH=LLH;
            initpars = [wG_hat kappa_hat];
        end
    end
end
fitpars = fminsearch(@(pars) mLLHfun(pars,x), initpars);

function mLLH = mLLHfun(pars,x)
wG = pars(1);
kappa = pars(2);
if wG<0 || wG>1 || kappa<0 || kappa>700
    mLLH = Inf;    
else
    mLLH = -sum(log(wG/2/pi + (1-wG)*1/2/pi/besseli(0,kappa)*exp(kappa*cos(x))));
end
