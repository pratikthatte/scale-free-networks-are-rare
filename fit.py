import numpy as np
import powerlaw as pl
import scipy.optimize as op
import scipy.special as sp
import time
def fit_power_law(x):
    fit = pl.Fit(x)
    ks = pl.power_law_ks_distance(x,fit.alpha,fit.xmin)
    return fit.xmin, fit.alpha, fit.n_tail,ks

def fit_power_law_and_compare_with_alternative(x):
    fit = pl.Fit(x)
    ks = pl.power_law_ks_distance(x,fit.alpha,fit.xmin)
    R_power_exponential, p_power_exponential = fit.distribution_compare('power_law','exponential')
    exponential_decision = decide(R_power_exponential,p_power_exponential)
    R_power_lognormal, p_power_lognormal = fit.distribution_compare('power_law','lognormal')
    lognormal_decision = decide(R_power_lognormal, p_power_lognormal)
    R_power_stretched, p_power_stretched = fit.distribution_compare('power_law','stretched_exponential')
    stretched_exponential_decision = decide(R_power_stretched, p_power_stretched)
    R_power_truncated, p_power_truncated = fit.distribution_compare('power_law','truncated_power_law')
    truncated_power_law_decision = decide_nested_power(R_power_truncated, p_power_truncated)
    return fit.xmin, fit.alpha, fit.n_tail,ks, exponential_decision, lognormal_decision, stretched_exponential_decision, truncated_power_law_decision
## Function to run goodness of fit test on the power law fit and calculate p-value. Taken from the fit.py class in this repository: https://github.com/adbroido/SFAnalysis
def power_law_pval(x, alpha, xmin, gof):
    eps = 0.01
    num_resamps = 1000
    bootstraps = np.zeros(num_resamps)
    n = len(x)
    xmax = np.max(x)
    tailinds = x>=xmin
    xtail = x[tailinds]
    xhead = x[~tailinds]
    ntail = len(xtail)
    nhead = len(xhead)
    ptail = float(ntail)/n
    mmax = 20*xmax
    const_tail = sp.zeta(alpha) - np.sum(np.arange(1,xmin)**(-alpha))
    pdf_tail = np.arange(xmin,mmax+1)**(-alpha)/const_tail
    pdf = np.zeros(mmax+1)
    pdf[xmin:] = pdf_tail
    del pdf_tail
    cdf = np.array( [ np.arange(mmax+1), np.cumsum(pdf) ] )
    cdf = np.concatenate( (cdf , np.array([[mmax+1,1]]).T) , axis = 1 )

    starttime = time.time()
    for resamp_ind in range(num_resamps):
        nnewhead = n
        while nnewhead >= n:
            nnewhead = np.sum(np.random.rand(n)>ptail)
        headinds = np.array([np.floor(nhead*np.random.rand(nnewhead))],dtype=int)
        newhead = xhead[headinds][0]
        nnewtail  = n-nnewhead

        # parametric bootstrap for the powerlaw tail
        rtail = np.sort(np.random.rand(nnewtail))
        newtail = np.zeros(nnewtail, dtype = int)
        indrtail = 0
        indnewtail = 0
        for xval in range(xmin, mmax+2):
            while (indrtail < len(rtail)) and (rtail[indrtail] <= cdf[1, xval]):
                indrtail += 1
            newtail[indnewtail:indrtail] = xval
            indnewtail = indrtail
            if indnewtail > nnewtail:
                break
        # combine into new sample
        newx = np.concatenate((newhead, newtail))
        if (newx == np.zeros_like(newx)).all():
            import pdb; pdb.set_trace()
        # fit this new sample
        [newxmin, newalpha, newntail, newgof] = fit_power_law(newx)
        current_p = np.sum(bootstraps[0:resamp_ind]>=gof)/(float(resamp_ind+1))
        bootstraps[resamp_ind] = newgof
        # if it's taking forever and we can end, do it
        if time.time() - starttime > 500:
            if resamp_ind > num_resamps/20.:
                if current_p<0.05 or current_p>0.5:
                    return current_p
    p = np.sum(bootstraps>=gof)/float(num_resamps)
    return p

def decide(R,p):
    if p <= 0.1:
        if R > 0:
            decision = 1
        elif R< 0:
            decision = -1
        else:
            decision = 0
    else:
        decision=0
    return decision

def decide_nested_power(R,p):
    if p <= 0.1 and R < 0 :
        decision = -1
    else:
        decision = 0
    return decision
