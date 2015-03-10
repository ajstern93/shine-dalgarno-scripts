from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import minimize
import random
from math import log

'''
Given a np.array of your data. Do:
results1, results2 = bimodal_test.get_gaussian_fits(data)
aic1_gauss aic2_gauss, oddsOf1_gauss = bimodal_test.compare_models(results1, results2)
results1, results2 = bimodal_test.get_poisson_fits(data)
aic1_poisson, aic2_poisson, oddsOf1_poisson = bimodal_test.compare_models(results1, results2)

If any model fits terribly, you will not get convergence on the optimal
parameters. I.e. aic1 and/or aic2 will be inf or nan. This just tells you
that your model is an abysmal fit.

'''

#######################################################################
#PDF functions, just really unpacking parameters. could probably be 
#written easier/shorter with *kwargs

def pdf_model_uni_poisson(x, params):
     lam = params
     return stats.poisson.pmf(x, lam)

def pdf_model_bi_poisson(x, params):
    lam1, lam2, pi_1 = params
    return pi_1*stats.poisson.pmf(x, lam1) + (1-pi_1)*stats.poisson.pmf(x, lam2)

def pdf_model_uni_gauss(x, params):
    mu, sig = params
    return stats.norm.pdf(x, mu, sig)

def pdf_model_bi_gauss(x, params):
    mu1, sig1, mu2, sig2, pi_1 = params
    return pi_1*stats.norm.pdf(x, mu1, sig1) + (1-pi_1)*stats.norm.pdf(x, mu2, sig2)

def log_likelihood_arbitrary(p, sample, function):
    '''
    likelihood of an arbitrary function
    may not work with exponential last time i tried due to zeroes
    so i have to investigate
    '''
    return -1*np.log(function(sample, p)).sum()


def get_gaussian_fits(data_1d):
    '''
    getting the best fitting gaussian and bimodal gaussian results
    '''
    p0_unimodal = np.array([np.mean(data_1d), np.std(data_1d)])
    results_unimodal = minimize(log_likelihood_arbitrary, x0=p0_unimodal, args=(data_1d,pdf_model_uni_gauss,), method='Nelder-Mead', options=dict(maxiter=10e3, maxfev=2e4))
    
    twenty_five_percentile = np.percentile(data_1d, 25)
    seventy_five_percentile = np.percentile(data_1d, 75)
    lower_half = [x for x in data_1d if x <= seventy_five_percentile]
    upper_half = [x for x in data_1d if x >= twenty_five_percentile]
    p0_bimodal = np.array([np.mean(lower_half), np.std(lower_half), np.mean(upper_half), np.std(upper_half), 0.5])
    results_bimodal = minimize(log_likelihood_arbitrary, x0=p0_bimodal, args=(data_1d,pdf_model_bi_gauss), method='Nelder-Mead', options=dict(maxiter=10e3, maxfev=2e4))
    
    return results_unimodal, results_bimodal

def get_poisson_fits(data_1d):
    '''
    getting the best fitting poisson and bimodal poisson results
    '''
    
    p0_unimodal = np.array([np.mean(data_1d)])
    results_unimodal = minimize(log_likelihood_arbitrary, x0=p0_unimodal, args=(data_1d,pdf_model_uni_poisson,), method='Nelder-Mead', options=dict(maxiter=10e3, maxfev=2e4))
    
    twenty_five_percentile = np.percentile(data_1d, 25)
    seventy_five_percentile = np.percentile(data_1d, 75)
    lower_half = [x for x in data_1d if x <= seventy_five_percentile]
    upper_half = [x for x in data_1d if x >= twenty_five_percentile]
    p0_bimodal = np.array([np.mean(lower_half), np.mean(upper_half), 0.5])
    results_bimodal = minimize(log_likelihood_arbitrary, x0=p0_bimodal, args=(data_1d,pdf_model_bi_poisson), method='Nelder-Mead', options=dict(maxiter=10e3, maxfev=2e4))
    
    return results_unimodal, results_bimodal

def compare_models(results1, results2):
    '''
    comparing models using AIC

    Third output is this: the probability that model 1 (aic1) is a better fit
    than model 2. So obviously keep track of what model you input as results 1
    and vice versa. 
    '''
    def get_AIC(likelihood, parameter_number):
            return (2*parameter_number) + 2*likelihood
    aic1 = get_AIC(results1.fun, len(results1.x))
    aic2 = get_AIC(results2.fun, len(results2.x))
    return aic1, aic2, np.exp((aic2-aic1)/2.)






