#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 13:28:16 2023

Dose-Response-Optim-Boot-NoIPy.py is a python3 script adapted from the code developed with Dr. Charles N. Haas 
earlier to conduct dose-response optimization in python as opposed to R as in CAMRA - Center for Advancing Microbial
Risk Assessment 

@author: Mark H. Weir Ph.D. - weir.95@osu.edu
"""
from typing import Sequence

"""Start by including the packages you will need"""
import json
import sys
from concurrent.futures import as_completed, ProcessPoolExecutor
from enum import Enum
from functools import partial
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as scipy
import scipy.special as sc
import scipy.stats
# These are for the ellipses
from matplotlib.patches import Ellipse
from scipy import optimize;


class Model(str, Enum):
    exp = 'exp'
    betapoisson = 'betapoisson'
    betabetapoisson = 'Betabetapoisson'
    exactbetapoisson = 'exactbetapoisson'


MODELS_TO_FIT = [
    Model.exp,
    Model.betapoisson,
    Model.betabetapoisson,
    Model.exactbetapoisson,
]


def get_cov_ellipse(cov, centre, nstd, **kwargs):
    """
    Function to use matplotlib Ellipse patch for the covariance matrix (cov)
    centred at centre and scaled by the factor nstd based on number of standard
    deviations.
    """

    # Find and sort eigenvalues and eigenvectors into descending order
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # The anti-clockwise angle to rotate our ellipse by 
    vx, vy = eigvecs[:, 0][0], eigvecs[:, 0][1]
    theta = np.arctan2(vy, vx)

    # Width and height of ellipse to draw
    width, height = 2 * nstd * np.sqrt(eigvals)
    return Ellipse(xy=centre, width=width, height=height,
                   angle=np.degrees(theta), **kwargs)


""" Establish the dose-response functions. Using the Exponential, beta Poisson, approximated
beta Poisson with alpha and N50 parameterization, and beta Poisson with alpha and beta parameterization.  """


def expDR(dose, k):
    prob = [1 - np.exp(-np.asarray(d) * k) for d in dose]
    return prob


# betapoissonDR - given vector of doses (d), parameter alpha, N50, return predicted probability of response
def betapoissonDR(dose, alpha, N50):
    prob = []
    z = -1 + 2 ** (1 / alpha)
    for d in dose:
        probn = (d / N50) * z
        probn = 1 - (1 + probn) ** (-alpha)
        prob.append(probn)
    return prob


def BetabetapoissonDR(dose, alpha, beta):
    prob = [1 - (1 + np.asarray(d) / beta) ** (-alpha) for d in dose]
    return prob


# exact beta poisson
def exactbetapoissonDR(dose, alpha, beta):
    prob = [1 - sc.hyp1f1(alpha, alpha + beta, -d) for d in dose]
    return prob


""" Wrapper function modifies the behavior of a the dose-response functions. This takes the 
 functions made and allows them to be used as arguments in other functions. """


def drwrapper(drmodel, drparameters, doses):
    if drmodel == 'exp':
        k = drparameters
        p = expDR(doses, k)
    elif drmodel == 'betapoisson':
        # param 1 = alpha, param 2=N50
        alpha, N50 = drparameters
        p = betapoissonDR(doses, alpha, N50)
    elif drmodel == 'Betabetapoisson':
        alpha, beta = drparameters
        p = BetabetapoissonDR(doses, alpha, beta)
    elif drmodel == 'exactbetapoisson':
        alpha, beta = drparameters
        p = exactbetapoissonDR(doses, alpha, beta)
    else:
        print("I Don't Recognize the Dose Response Model")
        sys.exit()
    return p


""" Deviance function --> This is where the wrapper is first used, note in the predprob
 the dose-response models are called as arguments to that function. Thus the wrapper 
 allows us to just have one deviance function defined and work for all of the dose-response
 models that we would want it to work for."""


def deviance(drparameters, drmodel, dataframe):
    # given model, call model function to get pred prob and return deviance
    # retrieve dose vector
    drparameters = np.exp(drparameters)  # remember inputs are log transformed
    doses = dataframe['dose']

    predprob = drwrapper(drmodel, drparameters, doses)
    Y = 0
    for i in dataframe.index:
        Ntot = dataframe.iloc[i]['subjects']
        Pos = dataframe.iloc[i]['positives']
        obsprob = Pos / Ntot
        predp = predprob[i]
        # evaluate Ypos --- check for pi_null=0
        Ypos = 0.0
        Yneg = 0.0
        if obsprob > 0:
            Ypos = Pos * np.log(predp / obsprob)
        # evaluate Yneg --- check for pi_null=1
        if obsprob < 1:
            Yneg = (Ntot - Pos) * np.log((1 - predp) / (1 - obsprob))
        Y = Y + Ypos + Yneg
    Y = -2 * Y
    return Y


# Operate Bootstrap and Plot Parameter Estimmates

def estimate(
        positives,
        dfboot: pd.DataFrame,
        modeltofit,
        deviance,
        fulldatabest,
        argstofunction,
        optimizationmethod,
        useoptions,
):
    dfboot['positives'] = positives
    # print(dfboot)
    argstofunction = (modeltofit, dfboot)

    results = optimize.minimize(
        deviance,
        fulldatabest,
        argstofunction,
        method=optimizationmethod,
        options=useoptions
    )
    return results


def run_dataset(
        df_dataset: pd.DataFrame = None,
        modeltofit: str = 'exactbetapoisson',
        optimizationmethod: str = 'Nelder-Mead',
        multicore: bool = False,
        dataset_name: str = '',
        dir_dataset: Path = Path.home()
):
    start_time = perf_counter()
    print(
        '============================================================================================================')
    print(f'Running {dataset_name}')
    print(
        '============================================================================================================')
    """ Set up what the pathogen and outcome is for the plots to be generated. Then import data and other 
     values that need to be included. The line df=... means that if you keep the dose-response optimization 
     script in the same directory as the data, then you will always have the correct pathname. """

    # source_path = Path(__file__).resolve()
    # source_dir = source_path.parent

    # pathogen='C. albicans'
    pathogen = df_dataset['Species']
    outcome = 'Mortality '

    if df_dataset is None:
        # df = pd.read_csv(source_dir/'field_CA30.csv')
        df = pd.read_csv('field_CA30.csv')
    else:
        df = df_dataset[['dose', 'positive_response', 'negative_response', 'n', 'obpos']]
        df.columns = ['dose', 'positives', 'negatives', 'subjects', 'obpos']
        df = df.reset_index()

    # modeltofit='exactbetapoisson'       #options are exp,  betapoisson (N50 approximation),
    # Betabetapoisson (beta approximation), exactbetapoisson
    # optimizationmethod = 'Nelder-Mead' #options include CG, Nelder-Mead
    bootstraptrials = 10000  # np.math.factorial(len(df)) #200
    guesses = [.1, .3]  # note this is logarithm of parameters.
    # Must be 1 param for exponential and 2 for others
    expguesses = -20
    print(df)

    """ Run the optimization """

    argstofunction = (modeltofit, df)
    useoptions = {'disp': True}

    if modeltofit == 'exp':
        results = optimize.minimize(deviance, expguesses, argstofunction,
                                    method=optimizationmethod, options=useoptions)
    else:
        results = optimize.minimize(deviance, guesses, argstofunction,
                                    method=optimizationmethod, options=useoptions)

    if modeltofit == 'exp':
        dof = len(df) - 1
    else:
        dof = len(df) - 2;

    gofstat = round(scipy.stats.chi2.ppf(0.95, dof), 4)
    gofPvalue = round(1 - scipy.stats.chi2.cdf(results.fun, dof), 4)
    bestfitChi = round(scipy.stats.chi2.ppf(0.95, 1), 4)

    """ Extract and organize the optimized parameters from the optimization , then output
    the estimated parameters dataframe."""

    if modeltofit == 'exp':
        MODEL = 'Exponential'
        LNparameters = results.x
        LNk = LNparameters[0]
        kOut = np.exp(LNk)
        ParametersOutput = pd.DataFrame({'Deviance': [results.fun], 'k': [kOut],
                                         'GOF Stat': [gofstat], 'GOF p-value': [gofPvalue],
                                         'Best Fit Stat': [bestfitChi]},
                                        columns=['Deviance', 'k', 'GOF Stat', 'GOF p-value', 'Best Fit Stat'])

    if modeltofit == 'betapoisson':
        MODEL = 'Approximated Beta Poisson'
        LNparameters = results.x
        LNalphaReport = LNparameters[0]
        LNN50Report = LNparameters[1]
        alphaReport = np.exp(LNalphaReport)
        N50Report = np.exp(LNN50Report)
        ParametersOutput = pd.DataFrame({'Deviance': [results.fun], 'alpha': [alphaReport], 'N50': [N50Report],
                                         'GOF Stat': [gofstat], 'GOF p-value': [gofPvalue],
                                         'Best Fit Stat': [bestfitChi]},
                                        columns=['Deviance', 'alpha', 'N50', 'GOF Stat', 'GOF p-value',
                                                 'Best Fit Stat'])

    if modeltofit == 'Betabetapoisson':
        MODEL = 'Approximated Beta Poisson'
        LNparameters = results.x
        LNalphaReport = LNparameters[0]
        LNbetaReport = LNparameters[1]
        alphaReport = np.exp(LNalphaReport)
        betaReport = np.exp(LNbetaReport)
        ParametersOutput = pd.DataFrame({'Deviance': [results.fun], 'alpha': [alphaReport], 'Approx beta': [betaReport],
                                         'GOF Stat': [gofstat], 'GOF p-value': [gofPvalue],
                                         'Best Fit Stat': [bestfitChi]},
                                        columns=['Deviance', 'alpha', 'Approx beta', 'GOF Stat', 'GOF p-value',
                                                 'Best Fit Stat'])

    if modeltofit == 'exactbetapoisson':
        MODEL = 'Exact Beta Poisson'
        LNparameters = results.x
        LNalphaReport = LNparameters[0]
        LNbetaReport = LNparameters[1]
        alphaReport = np.exp(LNalphaReport)
        betaReport = np.exp(LNbetaReport)
        ParametersOutput = pd.DataFrame({'Deviance': [results.fun], 'alpha': [alphaReport], 'beta': [betaReport],
                                         'GOF Stat': [gofstat], 'GOF p-value': [gofPvalue],
                                         'Best Fit Stat': [bestfitChi]},
                                        columns=['Deviance', 'alpha', 'beta', 'GOF Stat', 'GOF p-value',
                                                 'Best Fit Stat'])

    ParametersOutput.to_csv(dir_dataset / f'{modeltofit}-{outcome.strip().lower()}-parameters.csv', index=False)
    print(
        f'=====Results from Optimization for {modeltofit}=============================================================')
    print(ParametersOutput)

    """ Establish a predicted probability of response based on the model that was fit
     this will be useful for plotting later on. This is set up to have a finer granularity to 
     the dose data for the curves so that a spline or other assumed value function needs to be used
     for a smoothed line."""

    fulldatabest = results.x
    drparameters = np.exp(fulldatabest)
    # doses=df['dose']
    doses = np.geomspace((min(df['dose']) + 0.000001) / 3, max(df['dose']) * 3, num=50)
    bestpredprob = drwrapper(modeltofit, drparameters, doses)
    df['observed probability'] = df["positives"] / df['subjects']

    """ The plot the curve that you just made possible. No confidence intervals on this one. """

    plt.scatter(df['dose'], df['observed probability'], s=200)
    plt.plot(doses, bestpredprob, c='red', label=MODEL)
    plt.xscale('log')
    # plt.yscale('log')
    plt.grid(which='both', axis='both')
    # plt.title(pathogen+ outcome+ ' Dose-Response')
    # plt.legend()
    plt.xlabel('Dose')
    # plt.xlim(1000, pd.DataFrame.max(df['dose'])*2)
    plt.ylabel('Probability of Response')
    plt_path = dir_dataset / f'{modeltofit}-{outcome.strip().lower()}-dose-response.png'
    plt.savefig(plt_path)
    # plt.show()

    """ Bootstrap """

    # Set up Bootstrap
    dfboot = df.copy()  # duplicate dataframe and drop pred prob
    dfboot = dfboot.rename(columns={"positives": "original positives"})  # rename positives-->original positives

    dfboot['positives'] = 0
    rng = np.random.default_rng()
    bootstrap_parameters = np.empty((bootstraptrials, drparameters.size))

    i = 0;
    j = 0

    if multicore:
        with ProcessPoolExecutor() as executor:
            useoptions = {'disp': False}
            estimate_partial = partial(
                estimate,
                dfboot=dfboot,
                modeltofit=modeltofit,
                deviance=deviance,
                fulldatabest=fulldatabest,
                argstofunction=argstofunction,
                optimizationmethod=optimizationmethod,
                useoptions=useoptions
            )
            dfboot_subj_numpy = dfboot['subjects'].to_numpy()
            dfboot_obs_prob_numpy = dfboot['observed probability'].to_numpy()
            while i < bootstraptrials:
                # to draw binomial random vector use:  numpy.random.Generator.binomial
                futures = [
                    executor.submit(
                        estimate_partial,
                        positives=rng.binomial(dfboot_subj_numpy, dfboot_obs_prob_numpy)
                    )
                    for idx in range(bootstraptrials - i)
                ]

                for future in as_completed(futures):
                    results = future.result()
                    if results.success:
                        bootstrap_parameters[i, :] = results.x
                        i += 1
                        if i % 100 == 0:
                            print(f'{i=}, time={perf_counter() - start_time:.3f}s')
                    else:
                        j += 1
    else:
        while i < bootstraptrials:
            # to draw binomial random vector use:  numpy.random.Generator.binomial
            dfboot['positives'] = rng.binomial(dfboot['subjects'].to_numpy(),
                                               dfboot['observed probability'].to_numpy())
            # print(dfboot)
            argstofunction = (modeltofit, dfboot)
            useoptions = {'disp': False}
            results = optimize.minimize(deviance, fulldatabest, argstofunction,
                                        method=optimizationmethod, options=useoptions)
            if results.success:
                bootstrap_parameters[i, :] = results.x
                i += 1
                if i % 100 == 0:
                    print(f'{i=}, time={perf_counter() - start_time:.3f}s')
            else:
                j += 1

    print(f'dataset created in {perf_counter() - start_time:.3f}s')

    """ Outputs """

    """ Bootstrap Parameters CSVs """

    BootParms = np.exp(bootstrap_parameters)

    if modeltofit == 'exactbetapoisson':  # or 'Betabetapoisson':
        bp_output = pd.DataFrame(BootParms, columns=["alpha", "beta"])
    elif modeltofit == 'exp':
        bp_output = pd.DataFrame(BootParms, columns=["k"])
    elif modeltofit == 'betapoisson':
        bp_output = pd.DataFrame(BootParms, columns=['alpha', 'N50'])

    bp_ouput_name = f'{modeltofit}-{outcome.strip().lower()}-bootstrap-parameters.csv'
    bp_output_path = dir_dataset / bp_ouput_name
    bp_output.to_csv(bp_output_path, index=False)

    """ Plots """

    print(str(j) + " trials did not converge")
    base_plt_title = f'{pathogen} {outcome.strip()}, {MODEL}'
    plt_title = f'{base_plt_title} Bootstrapped Parameters'
    plt_path = dir_dataset / f'{modeltofit}-{outcome.strip().lower()}-bootstrap_params.png'
    # histogram (exponential k) or scatter plot (beta poisson)
    if modeltofit == 'exp':
        plt.hist(bootstrap_parameters, density=False)
        plt.xlabel('Ln(k)')
        plt.title(plt_title)
        plt.savefig(plt_path)
        plt.show()
        bootstrap_parameters = np.exp(bootstrap_parameters)
    elif modeltofit == 'betapoisson':
        # param 1 = alpha, param 2=N50
        bootstrap_parameters = np.exp(bootstrap_parameters)  # invert log transforms
        alpha = bootstrap_parameters[:, 0]
        N50 = bootstrap_parameters[:, 1]
        plt.scatter(alpha, N50)
        plt.title(plt_title)
        plt.xlabel('alpha')
        plt.ylabel('N50')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(which='both', axis='both')
        plt.savefig(plt_path)
        plt.show()
        # think about plt.hist2D
    elif modeltofit == 'Betabetapoisson':
        bootstrap_parameters = np.exp(bootstrap_parameters)  # invert log transforms
        alpha = bootstrap_parameters[:, 0]
        beta = bootstrap_parameters[:, 1]
        plt.scatter(alpha, beta)
        plt.title(plt_title)  # ', Bootstrap Trials = '+ str(bootstraptrials))
        plt.xlabel('alpha')
        plt.ylabel('beta')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(which='both', axis='both')
        plt.savefig(plt_path)
        plt.show()
    elif modeltofit == 'exactbetapoisson':
        bootstrap_parameters = np.exp(bootstrap_parameters)  # invert log transforms
        alpha = bootstrap_parameters[:, 0]
        beta = bootstrap_parameters[:, 1]
        plt.scatter(alpha, beta)
        plt.title(plt_title)  # ', Bootstrap Trials = '+ str(bootstraptrials))
        plt.xlabel('alpha')
        plt.ylabel('beta')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(which='both', axis='both')
        plt.savefig(plt_path)
        plt.show()
    else:
        print("I Don't Recognize the Dose Response Model")
        return

    """ Plot the dose-response curve with the bootstrapped confidence intervals """

    # recall - doses and best fits from before
    # so for each dose, calculate risk with all bootstrap params
    # then save 2.5 and 97.5 percentile
    bootstrapDR = pd.DataFrame()
    bootstrapDR["dose"] = doses
    bootstrapDR["bestpred"] = bestpredprob
    bootstrapDR["LCL"] = ""
    bootstrapDR["UCL"] = ""
    for i in bootstrapDR.index:
        d = [bootstrapDR.at[i, "dose"]]
        predprob = np.empty([bootstraptrials])
        for j in range(bootstraptrials):
            params = bootstrap_parameters[j, :]
            predprob[j] = drwrapper(modeltofit, params, d)[0]
        LCL, UCL = np.percentile(predprob, [2.5, 97.5])
        bootstrapDR.at[i, "LCL"] = LCL
        bootstrapDR.at[i, "UCL"] = UCL

    # do plot of dose response with bootstrap confidence regions
    plt.scatter(df['dose'], df['observed probability'], s=200)
    plt.plot(doses, bestpredprob, c='red', label='MLE of ' + MODEL)
    plt.plot(doses, bootstrapDR["LCL"], c='black', linestyle='dotted', label='')
    plt.plot(doses, bootstrapDR["UCL"], c='black', linestyle='dotted', label='')
    plt.xscale('log')
    # plt.yscale('log')
    plt.grid(which='both', axis='both')
    # plt.title(f'{base_plt_title} Model and Uncertainty')
    plt.xlabel('Dose')
    plt.ylabel('Probability of Response')
    # plt.legend()
    plt_path = dir_dataset / f'{modeltofit}-{outcome.strip().lower()}-model-and-uncertainty.png'
    plt.savefig(plt_path)
    # plt.show()
    return


def read_list_datasets_done(fp: Path) -> list[str]:
    datasets_done: list[str] = []
    if fp.exists():
        with open(fp, 'r') as f:
            datasets_done = json.load(f)
    return datasets_done


def write_list_datasets_done(datasets_done: list[str], fp: Path):
    with open(fp, 'w') as f:
        json.dump(datasets_done, f, indent=2)


def main(models_to_fit: Sequence[str] = MODELS_TO_FIT):
    dir_datasets = Path(r"C:\Users\steph\Documents\Maddy\dose-response")
    dir_datasets.mkdir(parents=True, exist_ok=True)
    file_datasets_todo = dir_datasets / 'candida_dr_data_master_all2.csv'
    file_datasets_done = dir_datasets / 'datasets_done.json'

    datasets_done = read_list_datasets_done(file_datasets_done)

    df_datasets_todo = pd.read_csv(file_datasets_todo)
    per_strain_list_df = df_datasets_todo[['Source', 'Species', 'Strain']].drop_duplicates().reset_index()
    for idx, row in per_strain_list_df.iterrows():

        df_dataset = df_datasets_todo[
            (df_datasets_todo['Source'] == row['Source']) &
            (df_datasets_todo['Species'] == row['Species']) &
            (df_datasets_todo['Strain'] == row['Strain'])
            ]
        dir_dataset = dir_datasets / 'autodata' / f"{row['Source']}_{row['Species']}_{row['Strain']}".replace(' ',
                                                                                                              '-').lower()
        dir_dataset.mkdir(parents=True, exist_ok=True)

        for model_to_fit in MODELS_TO_FIT:
            try:
                dataset_name = f"{row['Source']}_{row['Species']}_{row['Strain']}_{model_to_fit}".replace(' ',
                                                                                                          '-').lower()
                if dataset_name in datasets_done:
                    continue
                run_dataset(
                    df_dataset,
                    modeltofit=model_to_fit,
                    dataset_name=dataset_name,
                    multicore=False,
                    dir_dataset=dir_dataset
                )
                datasets_done.append(dataset_name)
                write_list_datasets_done(datasets_done, file_datasets_done)
            except Exception as e:
                print(f'!!!! Error running dataset, {dataset_name}')
                print(e)
    per_species_list_df = df_datasets_todo[['Source', 'Species']].drop_duplicates().reset_index()


if __name__ == '__main__':
    start_time = perf_counter()
    main()
    end_time = perf_counter()
    print(f'completed in {end_time - start_time:.3f}s')
# END OF CODE - the Navy can go to Hell
