"""Reproduce deterministic checks, independent Kepler simulations and figures.

Run from any directory: python validate_and_plot.py --n-per-e 400000
All observations generated here are synthetic; no survey data are analysed.
"""
from __future__ import annotations
import argparse
import json
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad, cumulative_simpson
from scipy.special import beta, fresnel
from scipy.stats import chi2, kstest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wide_binary_angles import (step_cdf, angle_pdf, powerlaw_angle_pdf, moment_observable,
    observable_bounds, estimate_moments, sample_kepler, odd_sine_moment,
    fractional_moment_observable, fractional_observable_bounds,
    eccentricity_square_observable, eccentricity_fourth_observable,
    powerlaw_mean_circularity, powerlaw_alpha_from_sin3, powerlaw_alpha_from_e2,
    powerlaw_fisher_information, abel_inverse_left_cdf,
    sample_kepler_physical, npmle_eccentricity)


def require(condition, message):
    """Explicit check that also runs under python -O (unlike assert)."""
    if not condition:
        raise RuntimeError("validation check failed: " + str(message))

ROOT = Path(__file__).resolve().parent
COLORS = ["#754568", "#B77642", "#427F72", "#667A98", "#9B5250"]
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10,
    "axes.spines.top":False, "axes.spines.right":False,
    "axes.labelcolor":"#343039", "text.color":"#343039",
    "axes.edgecolor":"#8F8792", "axes.prop_cycle":plt.cycler(color=COLORS),
    "figure.facecolor":"#FCFBF8", "axes.facecolor":"#FCFBF8",
    "savefig.facecolor":"#FCFBF8", "axes.titleweight":"semibold"})


def integral(e, func=lambda x:1.):
    peak = np.arccos(e)
    if e == 0:
        return 2*quad(lambda g: angle_pdf(g,e)*func(g),0,np.pi/2,
                      epsabs=1e-10, epsrel=1e-10, limit=300)[0]
    return 2*sum(quad(lambda g: angle_pdf(g,e)*func(g),a,b,
                     epsabs=1e-10, epsrel=1e-10,limit=300)[0]
                 for a,b in [(0,peak),(peak,np.pi/2)])


def _ideal_cdf_function(e, points=8001):
    """Vectorised P(Gamma<=x | e) on [0,pi/2]: exact quadrature on a fine grid
    (including the cusp) with linear interpolation; interpolation error is far
    below the Kolmogorov-Smirnov resolution at the sample sizes used."""
    peak=np.arccos(e) if e>0 else np.pi/2
    xs=np.unique(np.concatenate([np.linspace(0,np.pi/2,points),[peak]]))
    pieces=[quad(lambda g:angle_pdf(g,e),lo,hi,epsabs=1e-14,limit=50)[0]
            for lo,hi in zip(xs[:-1],xs[1:])]
    vals=2*np.concatenate([[0.],np.cumsum(pieces)])
    return lambda x:np.interp(x,xs,vals)


def deterministic_checks():
    rows=[]
    for e in [0.,1e-6,.1,.3,.6,.9,.99,.9999]:
        norm=integral(e)
        for k in [1,2,3,6,12]:
            actual=integral(e,lambda g:moment_observable(g,k))
            expected=(1-e*e)**(k-.5)
            rows.append(dict(e=e,k=k,normalisation=norm,
                             moment_integral=actual,exact=expected,
                             absolute_error=abs(actual-expected)))
    maxnorm=max(abs(r['normalisation']-1) for r in rows)
    maxmoment=max(r['absolute_error'] for r in rows)
    require(maxnorm<2e-8 and maxmoment<2e-8,(maxnorm,maxmoment))
    power=[]
    for alpha in [-.8,0.,1.,4.]:
        norm=2*quad(lambda g:powerlaw_angle_pdf(g,alpha),0,np.pi/2,
                    epsabs=1e-10)[0]
        for g in [.17,.61,1.2]:
            c=abs(np.cos(g))
            # Change variable u=e**(alpha+1) to remove an endpoint singularity.
            u0=c**(alpha+1)
            numeric=sum(quad(lambda u:angle_pdf(g,u**(1/(alpha+1))),a,b,
                             epsabs=1e-9,limit=300)[0]
                        for a,b in [(0,u0),(u0,1)])
            exact=powerlaw_angle_pdf(g,alpha)
            power.append(dict(alpha=alpha,gamma=g,normalisation=norm,
                              mixture_integral=numeric,closed_form=exact,
                              error=abs(numeric-exact)))
    require(max(r['error'] for r in power)<2e-7,'power-law mixture')
    # Near the centre and alpha=-1, sin(gamma)^2 rounds to one too early.
    # Evaluate the nonzero cosine directly in this sensitive limit.
    g=np.pi/2-1e-10;alpha_nc=-.999;alpha=alpha_nc
    reference=beta((alpha+1)/2,.5)/(2*np.pi)*(1-abs(np.cos(g))**(alpha+1))/(np.sin(g)**2)
    near_centre=powerlaw_angle_pdf(g,alpha);g_nc=g
    require(abs(near_centre/reference-1)<1e-12,'near-centre power law')
    endpoint=[]
    for e in [0.,.4,.9,.9999]:
        true=1/(2*np.pi*np.sqrt(1-e*e))
        actual=angle_pdf(1e-9,e)
        endpoint.append(dict(e=e,near_zero=actual,limit=true,
                             relative_error=abs(actual/true-1)))
    # End-to-end Abel inversion: recover F([0,r)) from the mixed angle density
    # alone (numerical derivative of Y, then the inverse formula). Atoms at
    # zero, .35 and .8, and continuous power laws, are both tested.
    atoms=[0.,.35,.8];weights=[.2,.5,.3]
    mixture=lambda g:sum(w*angle_pdf(g,e) for e,w in zip(atoms,weights))
    inversion=[]
    for r in [.1,.3,.5,.7,.9]:
        recovered=abel_inverse_left_cdf(mixture,r,atoms)
        expected=sum(w for e,w in zip(atoms,weights) if e<r)
        inversion.append(dict(law='atoms (0,.35,.8)',r=r,recovered_left_cdf=recovered,
                              exact=expected,error=abs(recovered-expected)))
    for alpha in [0.,1.,2.5]:
        for r in [.2,.6,.9]:
            recovered=abel_inverse_left_cdf(lambda g:powerlaw_angle_pdf(g,alpha),r)
            inversion.append(dict(law=f'power law alpha={alpha:g}',r=r,
                                  recovered_left_cdf=recovered,exact=r**(alpha+1),
                                  error=abs(recovered-r**(alpha+1))))
    require(max(r['error'] for r in inversion)<2e-3,'end-to-end Abel inversion')
    # Verify derivative of the combined kernel independently of the Abel beta integral.
    derivative=[]
    for e in [0.,.35,.8]:
        for c in [.12,.47,.92]:
            h=1e-5
            def y(t):return np.pi*(1-t*t)*angle_pdf(np.arccos(t),e)
            numerical=-(-y(c+2*h)+8*y(c+h)-8*y(c-h)+y(c-2*h))/(12*h)
            exact=1/np.sqrt(c*c-e*e) if c>e else 0.
            derivative.append(dict(e=e,c=c,numerical=numerical,exact=exact,
                                   error=abs(numerical-exact)))
    require(max(r['error'] for r in derivative)<2e-7,'kernel derivative')
    return dict(quadrature=rows,powerlaw=power,endpoints=endpoint,
                atomic_inversion=inversion,kernel_derivative=derivative,
                max_normalisation_error=maxnorm,max_moment_error=maxmoment,
                near_centre_powerlaw=dict(alpha=alpha_nc,gamma=g_nc,value=near_centre,reference=reference))


def physical_validation(n, rng):
    rows=[];histograms=[]
    edges=np.linspace(0,np.pi,25)
    for e in [0.,.3,.6,.9,.99]:
        state=sample_kepler(np.full(n,e),rng,return_state=True)
        g=state['gamma']
        momentrows=[]
        for k in range(1,7):
            v=moment_observable(g,k)
            target=(1-e*e)**(k-.5)
            se=v.std(ddof=1)/np.sqrt(n)
            momentrows.append(dict(k=k,estimate=float(v.mean()),exact=target,
                                    standard_error=float(se),z=float((v.mean()-target)/se)))
        masses=[]
        peak=np.arccos(e)
        for a,b in zip(edges[:-1],edges[1:]):
            points=[x for x in [peak,np.pi-peak] if a<x<b]
            masses.append(quad(lambda t:angle_pdf(t,e),a,b,points=points,
                               epsabs=1e-10,limit=200)[0])
        counts=np.histogram(g,edges)[0]
        expected=n*np.array(masses)
        stat=float(np.sum((counts-expected)**2/expected))
        # Binning-free Kolmogorov-Smirnov test against the exact CDF on [0,pi/2]
        # of the folded angle (the law is symmetric about pi/2).
        folded=np.minimum(g,np.pi-g)
        ks=kstest(folded,_ideal_cdf_function(e))
        rows.append(dict(e=e,n=n,kepler_solver_max_residual=state['kepler_residual'],
                         chi_square=stat,df=23,p_value=float(chi2.sf(stat,23)),
                         ks_statistic=float(ks.statistic),ks_p_value=float(ks.pvalue),
                         moments=momentrows))
        histograms.append((e,counts/n/np.diff(edges)))
    require(max(abs(z['z']) for r in rows for z in r['moments'])<5,'Monte Carlo moments')
    # Pre-registered gate: family-wise level 0.001 over ten distributional tests.
    require(min(min(r['p_value'],r['ks_p_value']) for r in rows)>1e-4,'distributional tests')
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    grid=np.linspace(.0005,np.pi-.0005,2400)
    for i,(e,hist) in enumerate(histograms):
        axes[0].plot(np.degrees(grid),angle_pdf(grid,e),color=COLORS[i],lw=1.7,label=f'e = {e:g}')
        axes[0].plot(np.degrees((edges[1:]+edges[:-1])/2),hist,'o',ms=3,
                     color=COLORS[i],alpha=.7)
    axes[0].set(xlabel='Projected angle Γ (degrees)',ylabel='Density per radian',
                title='Exact likelihood and independent Kepler simulations',ylim=(0,2.4),xlim=(0,180))
    axes[0].legend(frameon=False,ncol=2,fontsize=8)
    for a,col in zip([-.5,0,1,4],COLORS):
        axes[1].plot(np.degrees(grid),powerlaw_angle_pdf(grid,a),color=col,label=f'α = {a:g}')
    axes[1].set(xlabel='Projected angle Γ (degrees)',ylabel='Density per radian',
                title=r'Closed likelihood for $f(e) = (\alpha + 1)e^\alpha$',xlim=(0,180),ylim=(0,.82))
    axes[1].legend(frameon=False,fontsize=8)
    fig.savefig(ROOT/'figures'/'likelihood_validation.png',dpi=180)
    fig.savefig(ROOT/'figures'/'likelihood_validation.pdf')
    plt.close(fig)
    return rows


def moments_and_sensitivity(rng,n=400000):
    es=np.array([0.,.65,.95]); weights=np.array([.2,.5,.3])
    e=rng.choice(es,size=20000,p=weights)
    g=sample_kepler(e,rng)
    estimates=estimate_moments(g,6)
    truths=[float(np.dot(weights,(1-es*es)**(k-.5))) for k in range(1,7)]
    for row,true in zip(estimates['moments'],truths):row['population_truth']=true
    state=sample_kepler(np.full(n,.6),rng,return_state=True)
    g=state['gamma']
    selections=[('Ideal sampling',g),
        ('Inclination: |cos i| < 0.2',g[np.abs(state['cos_inclination'])<.2]),
        ('Inclination: |cos i| > 0.8',g[np.abs(state['cos_inclination'])>.8]),
        ('Projected separation > 1.1a',g[state['projected_separation']>1.1]),
        ('Uniform true anomaly',sample_kepler(np.full(n,.6),rng,phase='true'))]
    for deg in [10,30]:
        err=rng.normal(0,np.radians(deg),n)
        # Fold an additive signed direction error into the angle [0,pi].
        noisy=np.arccos(np.clip(np.cos(g+err),-1,1))
        selections.append((f'Angular error: σ = {deg}°',noisy))
    results=[]
    for label,angles in selections:
        v=moment_observable(angles,1)
        results.append(dict(scenario=label,n=len(v),true_circularity=.8,
                            estimate=float(v.mean()),bias=float(v.mean()-.8),
                            standard_error=float(v.std(ddof=1)/np.sqrt(len(v)))))
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained',
                          gridspec_kw={'width_ratios':[1,1.35]})
    kk=np.arange(1,7)
    means=[r['estimate'] for r in estimates['moments']]
    ses=[r['standard_error'] for r in estimates['moments']]
    axes[0].errorbar(2*kk-1,means,yerr=1.96*np.array(ses),fmt='o',color=COLORS[0],
                     capsize=3,label='Estimate ± 1.96 standard errors')
    axes[0].scatter(2*kk-1,truths,marker='x',color=COLORS[1],s=45,label='Population truth')
    axes[0].set(xlabel='Power of j = √(1 − e²)',ylabel='Population moment',
                title='Unbiased moments: 20,000 synthetic binaries',ylim=(-.04,1))
    axes[0].legend(frameon=False,fontsize=8)
    ypos=np.arange(len(results))
    axes[1].axvline(.8,color='#8B8690',lw=1,ls='--')
    axes[1].errorbar([r['estimate'] for r in results],ypos,
                     xerr=[1.96*r['standard_error'] for r in results],fmt='o',
                     color=COLORS[2],capsize=3)
    axes[1].set_yticks(ypos,[r['scenario'] for r in results],fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set(xlabel='Estimate of mean circularity',title='Sensitivity at fixed e = 0.6',xlim=(0,1.55))
    fig.savefig(ROOT/'figures'/'moments_and_sensitivity.png',dpi=180)
    fig.savefig(ROOT/'figures'/'moments_and_sensitivity.pdf')
    plt.close(fig)
    return dict(population=dict(eccentricities=es.tolist(),weights=weights.tolist()),
                benchmark=estimates,sensitivity=results)


def instability(delta=.75):
    cgrid=np.linspace(0,1,200001)
    gamma=np.linspace(0,np.pi,30001)
    c=np.abs(np.cos(gamma));s2=np.sin(gamma)**2
    es=np.linspace(0,1,2400)
    rows=[];curves={}
    for n in [1,2,4,8,16,32,64]:
        w=2*np.pi*n
        fs,fc=fresnel(np.sqrt(2*w/np.pi)*cgrid)
        a=np.sqrt(2*np.pi/w)*(np.cos(w*cgrid*cgrid)*fc+np.sin(w*cgrid*cgrid)*fs)
        integ=cumulative_simpson(a,x=cgrid,initial=0)
        yy=integ[-1]-integ
        yc=np.interp(c,cgrid,yy)
        with np.errstate(divide='ignore',invalid='ignore'):
            h=yc/(np.pi*s2)
        h[(gamma==0)|(gamma==np.pi)]=a[-1]/(2*np.pi)
        tv=float(delta*np.trapezoid(np.abs(h),gamma))
        pplus=1/np.pi+delta*h;pminus=1/np.pi-delta*h
        kl=float(np.trapezoid(pplus*np.log(pplus/pminus),gamma))
        rows.append(dict(n=n,latent_total_variation=2*delta/np.pi,
                         angle_total_variation=tv,angle_kl_plus_minus=kl,
                         plus_normalisation=float(np.trapezoid(pplus,gamma))))
        curves[n]=(pplus,pminus)
    fig,axes=plt.subplots(2,2,figsize=(10.2,6.2),layout='constrained')
    for rr,n in enumerate([2,32]):
        for sign,col in [(1,COLORS[0]),(-1,COLORS[1])]:
            axes[rr,0].plot(es,2*es*(1+sign*delta*np.cos(2*np.pi*n*es*es)),color=col,lw=1.3)
        axes[rr,0].set(xlabel='Eccentricity e',ylabel='Population density',title=f'n = {n}; latent TV = {2*delta/np.pi:.3f}',xlim=(0,1))
        for curve,col in zip(curves[n],COLORS):
            axes[rr,1].plot(np.degrees(gamma),curve,color=col,lw=1.3)
        row=next(r for r in rows if r['n']==n)
        axes[rr,1].set(xlabel='Projected angle Γ (degrees)',ylabel='Density per radian',
                       title=f'Observed TV = {row["angle_total_variation"]:.4f}',xlim=(0,180),ylim=(.22,.42))
    fig.savefig(ROOT/'figures'/'inverse_instability.png',dpi=180)
    fig.savefig(ROOT/'figures'/'inverse_instability.pdf')
    plt.close(fig)
    return dict(delta=delta,method='Fresnel evaluation of Abel transform; Simpson integration on 200001 c points',rows=rows)


def _ideal_expectation(func, e):
    """E[func(Gamma)|e] by deterministic quadrature, split at the kernel cusp."""
    peak=np.arccos(e) if e>0 else np.pi/2
    return 2*sum(quad(lambda g:angle_pdf(g,e)*func(g),a,b,epsabs=1e-12,
                      epsrel=1e-11,limit=400)[0] for a,b in [(0,peak),(peak,np.pi/2)])


def _powerlaw_expectation(func, alpha):
    return 2*quad(lambda g:powerlaw_angle_pdf(g,alpha)*func(g),0,np.pi/2,
                  epsabs=1e-13,epsrel=1e-11,limit=400)[0]


def extensions(rng, n=400000, replicates=400, n_rep=5000):
    """Checks for the fractional-moment, Fourier and power-law-index results."""
    frac=[]
    for e in [0.,.3,.6,.9,.99]:
        j=np.sqrt(1-e*e)
        for s in [.25,.5,1.,2.,3.5,6.]:
            v=_ideal_expectation(lambda g:fractional_moment_observable(g,s),e)
            frac.append(dict(e=e,s=s,quadrature=v,exact=j**s,error=abs(v-j**s)))
        for name,obs,target in [('e^2',eccentricity_square_observable,e**2),
                                ('e^4',eccentricity_fourth_observable,e**4)]:
            v=_ideal_expectation(obs,e)
            frac.append(dict(e=e,s=name,quadrature=v,exact=target,error=abs(v-target)))
    require(max(r['error'] for r in frac)<1e-8,'fractional and Fourier moment identities')
    # Bounds: dense-grid check of the stated ranges.
    grid=np.linspace(0,np.pi,200001)
    bounds=[]
    for s in [0.,.5,1.,2.,3.5,6.]:
        lo,hi=fractional_observable_bounds(s)
        v=fractional_moment_observable(grid,s)
        bounds.append(dict(s=s,stated=[lo,hi],grid=[float(v.min()),float(v.max())]))
        require(v.min()>=lo-1e-12 and v.max()<=hi+1e-12 and v.min()-lo<1e-6 and hi-v.max()<1e-9,f'bounds s={s}')
    for obs,lo,hi in [(eccentricity_square_observable,-1.,1.25),
                      (eccentricity_fourth_observable,1/3-2*np.sqrt(2)/3,1/3+2*np.sqrt(2)/3)]:
        v=obs(grid);require(v.min()>=lo-1e-12 and v.max()<=hi+1e-12 and v.min()-lo<1e-6 and hi-v.max()<1e-6,'Fourier observable range')
    # Independent Kepler simulation of the new observables.
    mc=[]
    for e in [0.,.6,.95]:
        g=sample_kepler(np.full(n,e),rng)
        for name,vals,target in [('j^0.5',fractional_moment_observable(g,.5),(1-e*e)**.25),
                                 ('j^2',fractional_moment_observable(g,2.),1-e*e),
                                 ('e^2',eccentricity_square_observable(g),e*e),
                                 ('e^4',eccentricity_fourth_observable(g),e**4)]:
            se=float(vals.std(ddof=1)/np.sqrt(n))
            mc.append(dict(e=e,observable=name,estimate=float(vals.mean()),exact=target,
                           standard_error=se,z=float((vals.mean()-target)/se)))
    require(max(abs(r['z']) for r in mc)<6,'Monte Carlo fractional/Fourier moments')
    # Fisher information and moment-estimator efficiency for the power-law index.
    W1=lambda g:moment_observable(g,1)
    fisher=[]
    for alpha in [-.5,0.,.5,1.,1.3,2.,3.]:
        info=powerlaw_fisher_information(alpha)
        h=1e-5
        dmj=(powerlaw_mean_circularity(alpha+h)-powerlaw_mean_circularity(alpha-h))/(2*h)
        vj=_powerlaw_expectation(lambda g:W1(g)**2,alpha)-powerlaw_mean_circularity(alpha)**2
        me2=(alpha+1)/(alpha+3);dme2=2/(alpha+3)**2
        ve2=_powerlaw_expectation(lambda g:(1-eccentricity_square_observable(g))**2,alpha)-(1-me2)**2
        fisher.append(dict(alpha=alpha,fisher_per_angle=info,fisher_per_known_e=1/(alpha+1)**2,
            information_ratio=info*(alpha+1)**2,mle_sd_times_sqrtN=1/np.sqrt(info),
            sin3_estimator_sd_times_sqrtN=float(np.sqrt(vj)/abs(dmj)),
            e2_estimator_sd_times_sqrtN=float(np.sqrt(ve2)/dme2),
            sin3_efficiency=float((1/info)/(vj/dmj**2)),e2_efficiency=float((1/info)/(ve2/dme2**2))))
    require(all(0<r['sin3_efficiency']<=1+1e-9 and 0<r['e2_efficiency']<=1+1e-9 for r in fisher),'Cramer-Rao consistency')
    # Efficiency of the sin^3 estimator on a 501-point grid over [0.5, 3].
    grid_eff=[]
    for alpha in np.linspace(.5,3,501):
        h=1e-5
        dmj=(powerlaw_mean_circularity(alpha+h)-powerlaw_mean_circularity(alpha-h))/(2*h)
        vj=_powerlaw_expectation(lambda g:W1(g)**2,alpha)-powerlaw_mean_circularity(alpha)**2
        grid_eff.append((1/powerlaw_fisher_information(alpha))/(vj/dmj**2))
    grid_min=dict(alpha_range=[.5,3.],points=501,min_sin3_efficiency=float(min(grid_eff)),
                  argmin_alpha=float(np.linspace(.5,3,501)[int(np.argmin(grid_eff))]))
    require(grid_min['min_sin3_efficiency']>.95,'sin^3 efficiency grid')
    # Replicate experiment at the thermal law: sampling spread of the two estimators.
    est_s3=[];est_e2=[]
    for _ in range(replicates):
        e=np.sqrt(rng.uniform(size=n_rep))
        g=sample_kepler(e,rng)
        est_s3.append(powerlaw_alpha_from_sin3(g));est_e2.append(powerlaw_alpha_from_e2(g))
    est_s3=np.array(est_s3);est_e2=np.array(est_e2)
    thermal=next(r for r in fisher if r['alpha']==1.)
    replicate=dict(alpha=1.,n_per_sample=n_rep,replicates=replicates,
        sin3_mean=float(np.nanmean(est_s3)),sin3_sd=float(np.nanstd(est_s3,ddof=1)),
        e2_mean=float(np.nanmean(est_e2)),e2_sd=float(np.nanstd(est_e2,ddof=1)),
        predicted_sin3_sd=thermal['sin3_estimator_sd_times_sqrtN']/np.sqrt(n_rep),
        predicted_e2_sd=thermal['e2_estimator_sd_times_sqrtN']/np.sqrt(n_rep),
        predicted_cramer_rao_sd=thermal['mle_sd_times_sqrtN']/np.sqrt(n_rep),
        undefined_sin3=int(np.isnan(est_s3).sum()),undefined_e2=int(np.isnan(est_e2).sum()))
    require(abs(replicate['sin3_sd']/replicate['predicted_sin3_sd']-1)<.15,'replicate spread (sin^3)')
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    al=[r['alpha'] for r in fisher]
    axes[0].plot(al,[r['information_ratio'] for r in fisher],'o-',color=COLORS[0],label='Angle / known eccentricity')
    axes[0].plot(al,[r['sin3_efficiency'] for r in fisher],'s-',color=COLORS[2],label='sin³ moment estimator efficiency')
    axes[0].plot(al,[r['e2_efficiency'] for r in fisher],'^-',color=COLORS[1],label='⟨e²⟩ estimator efficiency')
    axes[0].set(xlabel='Power-law index α in f(e) = (α + 1)e^α',ylabel='Ratio',ylim=(0,1.05),
                title='Information in one projected angle')
    axes[0].legend(frameon=False,fontsize=8)
    bins=np.linspace(.5,1.5,41)
    axes[1].hist(est_s3,bins=bins,color=COLORS[2],alpha=.65,label='From mean sin³Γ')
    axes[1].hist(est_e2,bins=bins,color=COLORS[1],alpha=.45,label='From ⟨e²⟩ observable')
    axes[1].axvline(1,color='#8B8690',lw=1,ls='--')
    axes[1].set(xlabel='Estimated α',ylabel='Replicate samples',
                title=f'{replicates} synthetic thermal samples of {n_rep:,} binaries')
    axes[1].legend(frameon=False,fontsize=8)
    fig.savefig(ROOT/'figures'/'powerlaw_index_information.png',dpi=180)
    fig.savefig(ROOT/'figures'/'powerlaw_index_information.pdf')
    plt.close(fig)
    return dict(identities=frac,bounds=bounds,kepler_monte_carlo=mc,fisher=fisher,replicate=replicate,
                sin3_efficiency_grid=grid_min,
                max_identity_error=max(r['error'] for r in frac),
                max_mc_z=max(abs(r['z']) for r in mc))


def _semimajor_axes(rng, n, beta, amin=10., amax=1e6):
    """Draw a with density proportional to a**(-beta) on [amin, amax]."""
    u=rng.uniform(size=n)
    if beta==1:
        return amin*(amax/amin)**u
    k=1-beta
    return (amin**k+u*(amax**k-amin**k))**(1/k)


def selection_and_estimation(rng, n=2000000):
    """Phase controls, azimuth-blind thermal invariance, separation binning
    under power-law semimajor axes, and a grid NPMLE demonstration."""
    edges=np.linspace(0,np.pi,25)
    def chi_vs_kernel(g,e):
        peak=np.arccos(e)
        masses=[quad(lambda t:angle_pdf(t,e),a,b,points=[x for x in [peak,np.pi-peak] if a<x<b] or None,
                     epsabs=1e-11,limit=200)[0] for a,b in zip(edges[:-1],edges[1:])]
        counts=np.histogram(g,edges)[0];expected=len(g)*np.array(masses)
        stat=float(np.sum((counts-expected)**2/expected));return stat,float(chi2.sf(stat,23))
    controls=[]
    for e in [.6,.9]:
        ge=sample_kepler_physical(np.full(400000,e),np.ones(400000),rng,phase='eccentric')['gamma']
        gf=sample_kepler(np.full(400000,e),rng,phase='true')
        se=chi_vs_kernel(ge,e);sf=chi_vs_kernel(gf,e)
        controls.append(dict(e=e,uniform_eccentric_anomaly=dict(chi_square=se[0],p_value=se[1]),
                             uniform_true_anomaly=dict(chi_square=sf[0],p_value=sf[1])))
    require(all(c['uniform_eccentric_anomaly']['p_value']>1e-4 for c in controls),'uniform-E control must agree')
    require(all(c['uniform_true_anomaly']['p_value']<1e-10 for c in controls),'uniform-f control must be rejected')
    # Thermal invariance under azimuth-blind selections (log-flat a).
    e=np.sqrt(rng.uniform(size=n));a=_semimajor_axes(rng,n,1.)
    st=sample_kepler_physical(e,a,rng)
    s=st['projected_separation'];v=st['projected_speed']
    inbin=(s>1e3)&(s<10**3.5)
    fast=inbin&(v*np.sqrt(s)>np.median((v*np.sqrt(s))[inbin]))
    flat=lambda g:float(chi2.sf(np.sum((np.histogram(g,edges)[0]-len(g)/24)**2/(len(g)/24)),23))
    thermal=dict(n_bin=int(inbin.sum()),p_flat_bin=flat(st['gamma'][inbin]),
                 n_bin_fast=int(fast.sum()),p_flat_bin_fast=flat(st['gamma'][fast]))
    require(thermal['p_flat_bin']>1e-4 and thermal['p_flat_bin_fast']>1e-4,'thermal invariance')
    # Separation binning with dN/da ~ a^-beta; compare W_1 with the parent mean
    # circularity and with the exact rho^(beta-1) reweighting prediction.
    binning=[]
    for alpha in [0.,1.3]:
        for beta in [1.,1.6,2.]:
            e=rng.uniform(size=n)**(1/(alpha+1));a=_semimajor_axes(rng,n,beta)
            st=sample_kepler_physical(e,a,rng);s=st['projected_separation'];g=st['gamma']
            inbin=(s>1e3)&(s<10**3.5);rho=s/a;W=moment_observable(g,1)
            weight=rho**(beta-1)
            binning.append(dict(alpha=alpha,beta=beta,n_bin=int(inbin.sum()),
                parent_mean_circularity=powerlaw_mean_circularity(alpha),
                selected_mean_circularity=float(np.sqrt(1-e[inbin]**2).mean()),
                w1_estimate=float(W[inbin].mean()),w1_standard_error=float(W[inbin].std(ddof=1)/np.sqrt(inbin.sum())),
                reweighting_prediction=float(np.sum(W*weight)/np.sum(weight)),
                alpha_hat_sin3=powerlaw_alpha_from_sin3(g[inbin])))
    opik=[r for r in binning if r['beta']==1.]
    require(all(abs(r['w1_estimate']-r['parent_mean_circularity'])<5*r['w1_standard_error'] for r in opik),'Opik invariance')
    # Grid NPMLE demonstration on a mixture with a uniform bump, thermal part and an atom.
    def draw(m):
        u=rng.uniform(size=m);x=np.empty(m)
        b1=u<.3;x[b1]=rng.uniform(.1,.4,b1.sum())
        b2=(u>=.3)&(u<.8);x[b2]=np.sqrt(rng.uniform(size=b2.sum()))
        x[u>=.8]=.9;return x
    Ftrue=lambda r:.3*np.clip((r-.1)/.3,0,1)+.5*r*r+.2*(r>=.9)
    rs=np.array([.2,.3,.5,.7,.85,.95]);npmle=[]
    for m in [2000,20000,100000]:
        x=draw(m);g=sample_kepler(x,rng)
        fit=npmle_eccentricity(g);grid,wt=fit['grid'],fit['weights']
        est=step_cdf(grid,wt,rs)
        npmle.append(dict(n=m,r=rs.tolist(),cdf_estimate=est.tolist(),cdf_true=Ftrue(rs).tolist(),
                          max_abs_cdf_error=float(np.max(np.abs(est-Ftrue(rs)))),
                          mean_circularity_estimate=float(np.dot(wt,np.sqrt(1-grid**2))),
                          sample_mean_circularity_latent=float(np.mean(np.sqrt(1-x*x))),
                          w1_mean_circularity=float(np.mean(moment_observable(g,1))),
                          loglik_per_obs=fit['loglik_per_obs'],iterations=fit['iterations'],
                          optimality_residual=fit['optimality_residual'],
                          capped_kernel_entries=fit['capped_kernel_entries'],
                          weights=wt.tolist()))
    require(npmle[-1]['max_abs_cdf_error']<.06,'NPMLE demonstration')
    return dict(phase_controls=controls,thermal_invariance=thermal,separation_binning=binning,npmle=npmle)


def main():
    global ROOT
    parser=argparse.ArgumentParser();parser.add_argument('--n-per-e',type=int,default=400000)
    parser.add_argument('--out',default=None,help='output directory (default: alongside this script)')
    args=parser.parse_args()
    if args.out:
        ROOT=Path(args.out).resolve()
    (ROOT/'results').mkdir(parents=True,exist_ok=True);(ROOT/'figures').mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();rng=np.random.default_rng(20261028)
    out=dict(seed=20261028,python=platform.python_version(),numpy=np.__version__,
             scipy=scipy.__version__,synthetic_data_only=True)
    out['deterministic']=deterministic_checks()
    out['physical_validation']=physical_validation(args.n_per_e,rng)
    out['moments_and_sensitivity']=moments_and_sensitivity(rng,args.n_per_e)
    out['inverse_instability']=instability()
    out['extensions']=extensions(rng,args.n_per_e)
    out['selection_and_estimation']=selection_and_estimation(rng)
    out['elapsed_seconds']=time.perf_counter()-start
    (ROOT/'results'/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
    summary=dict(max_normalisation_error=out['deterministic']['max_normalisation_error'],
        max_moment_error=out['deterministic']['max_moment_error'],
        maximum_mc_moment_z=max(abs(z['z']) for r in out['physical_validation'] for z in r['moments']),
        histogram_checks=[{k:r[k] for k in ['e','n','chi_square','df','p_value','ks_statistic','ks_p_value']} for r in out['physical_validation']],
        sensitivity=out['moments_and_sensitivity']['sensitivity'],
        instability=out['inverse_instability']['rows'],
        abel_inversion_max_error=max(r['error'] for r in out['deterministic']['atomic_inversion']),
        extension_max_identity_error=out['extensions']['max_identity_error'],
        extension_max_mc_z=out['extensions']['max_mc_z'],
        powerlaw_fisher=out['extensions']['fisher'],
        sin3_efficiency_grid=out['extensions']['sin3_efficiency_grid'],
        thermal_replicate=out['extensions']['replicate'],
        selection_and_estimation=out['selection_and_estimation'],
        elapsed_seconds=out['elapsed_seconds'])
    summary['pooled_chi_square']=sum(r['chi_square'] for r in out['physical_validation'])
    summary['pooled_df']=23*len(out['physical_validation'])
    summary['pooled_p_value']=float(chi2.sf(summary['pooled_chi_square'],summary['pooled_df']))
    (ROOT/'results'/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
