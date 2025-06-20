import marimo

__generated_with = "0.7.5"
app = marimo.App(width="full")


@app.cell
def __():
    ## import all packages here
    import marimo as mo

    # data analysis packages
    import numpy as np
    import pandas as pd
    import polars as pl
    import random
    from collections import Counter
    from functools import reduce

    # stats packages
    from scipy.stats import poisson
    from scipy.stats import norm
    from scipy import stats
    import scipy.optimize as optim
    from scipy.optimize import curve_fit

    # plotting packages
    import altair as alt
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.pyplot import cm

    # ode package
    from scipy.integrate import odeint, solve_ivp # DE solver
    return (
        Counter,
        alt,
        cm,
        curve_fit,
        mo,
        norm,
        np,
        odeint,
        optim,
        pd,
        pl,
        plt,
        poisson,
        random,
        reduce,
        sns,
        solve_ivp,
        stats,
    )


@app.cell
def __():
    # The maximum allowed step size for the ODE solver.
    # Decreasing it won't decrease behavior, increasing it might increase numerical error.
    # 0.1 is a good reference value.
    t_step = 0.1
    return t_step,


@app.cell
def __():
    figureWidth = 150
    figureHeight = 150
    return figureHeight, figureWidth


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Here I made a matrix of growth rate given both the atb concentrations and the TCN.  
        The idea is that:  
        A. Without atb selections, i.e., [atb] = 0:  
            i. the cells that do not contain any transposon, i.e., TCN = 0, has the highest growth rate  
            ii. increasing TCN linearly decreases growth rate by a constant factor  
        B. With atb selections, i.e., [atb] > 0:  
            i. TCN = 0 -> cell death  
            ii. increasing TCN increases MIC  
                We implement this increase using hill equation by increasing the [atb] needed to reduce the growth rate by half. We term this as HMIC.
                $u_i = u_{max}(1-\beta i)\frac{[(1+\alpha i)HMIC]^2}{A^2+(1+\alpha i)HMIC]^2}$

        ## The growth rates of subpopulations determine the growth advantages, which determines the outcome of the TCN. It's the most important function in this model.
        """
    )
    return


@app.cell
def __():
    ATB_CONCS = [0, 0.5, 1, 2, 3]
    return ATB_CONCS,


@app.cell
def __(np):
    def calculate_growth_rate_matrix(antibiotic_concs, transposon_copies, baseline_mu, transposon_linear_decrease_rate = 0.003, MIC_increase_factor = 0.03):
        ### Use 0.003 & 0.03 for main figures; 0.005 & 0.05/0.08 for supplementary fig.10
        """
        Calculate a matrix of growth rates based on antibiotic concentrations and transposon copy numbers.

        Parameters:
        antibiotic_concs (list or array): List of antibiotic concentrations
        transposon_copies (list or array): List of transposon copy numbers

        Returns:
        numpy.ndarray: 2D array (matrix) of growth rates
        """
        base_growth_rate = baseline_mu  # Max baseline growth rate with plasmids, without transposons

        # Constants for Hill equation (resistance effect)
        hill_coeff = 2.0  # Hill coefficient (steepness of the curve)
        MIC = 0.001

        # Convert inputs to numpy arrays
        antibiotic_concs = np.array(antibiotic_concs)
        transposon_copies = np.array(transposon_copies)

        # Create meshgrid for vectorized calculation
        antibiotic_mesh, transposon_mesh = np.meshgrid(antibiotic_concs, transposon_copies)

        # If transposon is resistant, more transposon means higher MIC, less effective drug effect, linear increase from the baseline MIC
        transposon_MIC_mesh = MIC+(MIC_increase_factor*transposon_mesh)
        resistance_effect = transposon_MIC_mesh**hill_coeff / (antibiotic_mesh**hill_coeff + transposon_MIC_mesh**hill_coeff)

        # Calculate linear decrease effect (less pronounced)
        transposon_linear_decrease = 1 - transposon_linear_decrease_rate * transposon_mesh

        # Calculate final growth rate matrix
        growth_rate_matrix = base_growth_rate * resistance_effect * transposon_linear_decrease

        # Set the cells without any transposon and with no selection back to baseline growth rate
        growth_rate_matrix[0,0] = base_growth_rate

        return growth_rate_matrix
    return calculate_growth_rate_matrix,


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Here I made a function to call the calculate_growth_rate_matrix with defined max growth rate and PCN.
        It's useful since max TCN = PCN + 1.  
        This function takes care of this additional layer of logic, since we need to vary the PCN in simulation, given the plasmid amplification property we are looking at.  
        It's called growthRateFluctuating since the experimental fluctuating atb treatment fluctuates in between 0 and max/very high [atb] (or lower [atb] if needed). This matrix includes all possible growth rates.
        """
    )
    return


@app.cell
def __(ATB_CONCS, calculate_growth_rate_matrix, np):
    ### It takes 2 inputs:

    ### plasmidNum: the number of plasmids;
    ###          If there are n plasmids, then there are n+1 maximum possible subpopulations,
    ###          each subpopulation has a different number of transposons on plasmids, from 0 to n

    ### It returns growth rates for all subpopulations under two conditions
    ### i.e., selection or no selection; selection A or selection B

    def growthRateFluctuating(plasmidNum, baselineMu):
        # The matrix to be returned that
        muMatrix = []
        # A list of all possible transposon numbers
        transposonNumList = np.arange(0, plasmidNum+1, 1)  
        muMatrix = calculate_growth_rate_matrix(ATB_CONCS, transposonNumList, baselineMu) 

        return muMatrix
    return growthRateFluctuating,


@app.cell
def __():
    return


@app.cell
def __(growthRateFluctuating, np):
    ### Testing here
    popSizeList = np.arange(10, 101, 5)
    muBaselineList = 1*np.ones(len(popSizeList))

    plasmidN01, muBaseline01 = popSizeList[0], muBaselineList[0]
    plasmidN02, muBaseline02 = popSizeList[9], muBaselineList[9]
    plasmidN03, muBaseline03 = popSizeList[-1], muBaselineList[-1]

    # Run
    growthRateMatrix_test01 = growthRateFluctuating(plasmidN01, muBaseline01)
    growthRateMatrix_test02 = growthRateFluctuating(plasmidN02, muBaseline02)
    growthRateMatrix_test03 = growthRateFluctuating(plasmidN03, muBaseline03)

    growthRateMatrix_test = [growthRateMatrix_test03]
    plasmidList = [plasmidN03]
    return (
        growthRateMatrix_test,
        growthRateMatrix_test01,
        growthRateMatrix_test02,
        growthRateMatrix_test03,
        muBaseline01,
        muBaseline02,
        muBaseline03,
        muBaselineList,
        plasmidList,
        plasmidN01,
        plasmidN02,
        plasmidN03,
        popSizeList,
    )


@app.cell
def __(growthRateMatrix_test, np, plasmidList, plt):
    ## Sanity check: Plot the current growth rate
    # Get the plasma colormap
    cmap = plt.get_cmap('Greens_r')
    color1 = cmap(0.3)  # 50% of the way through the colormap
    color = [color1]

    ### Figure for one growth curves under selection
    fig = plt.figure(figsize=(2,2))
    for g in range(len(growthRateMatrix_test)): 
        # print(g)
        c = color[g]
        # print(c)
        PCN = plasmidList[g]
        LineLabel = 'PCN = ' + str(PCN)
        growthRateMatrix = growthRateMatrix_test[g]
        plt.plot(range(len(growthRateMatrix.T[0])), growthRateMatrix.T[0], linewidth = 3.0, color=c, label=LineLabel)
        # plt.title('growth rate under no selection', fontsize = 20)
        plt.ylim([-0.1, 1.1])
        plt.xlim([-10, 110])
        # plt.ylim([-0.1, 1.3])
        # plt.xticks(fontsize=40)
        plt.xticks([0, 100], fontsize=20)
        # plt.yticks(fontsize=40)
        plt.yticks(np.arange(0, 1.1, 0.5), fontsize=20)
        # Remove tick lines while keeping labels
        plt.tick_params(axis='both', which='both', length=0)
    plt.show()
    # fig.savefig("Fig2DLeft.pdf", bbox_inches="tight")

    ### Figure for both growth curves with the curves of mu_noS dimmed
    fig = plt.figure(figsize=(2,2))
    for g in range(len(growthRateMatrix_test)): 
        # print(g)
        c = color[g]
        # print(c)
        PCN = plasmidList[g]
        LineLabel = 'PCN = ' + str(PCN)
        growthRateMatrix = growthRateMatrix_test[g]
        plt.plot(range(len(growthRateMatrix.T[0])), growthRateMatrix.T[0], linewidth = 3.0, color=c, alpha = 0.15) #, label=LineLabel
        # plt.title('growth rate under no selection', fontsize = 20)
        plt.ylim([-0.1, 1.1])
        plt.xlim([-10, 110])
        # plt.ylim([-0.1, 1.3])
        # plt.xticks(fontsize=40)
        plt.xticks([0, 100], fontsize=20)
        # plt.yticks(fontsize=40)
        plt.yticks(np.arange(0, 1.1, 0.5), fontsize=20)
        # Remove tick lines while keeping labels
        plt.tick_params(axis='both', which='both', length=0)

    # fig = plt.figure()
    for g in range(len(growthRateMatrix_test)):  
        c = color[g]
        # print(c)
        PCN = plasmidList[g]
        LineLabel = 'PCN = ' + str(PCN)
        growthRateMatrix = growthRateMatrix_test[g]
        plt.plot(range(len(growthRateMatrix.T[-1])), growthRateMatrix.T[-1], linewidth = 3.0, color=c, label=LineLabel)
        # plt.title('growth rate under selection', fontsize = 20)
        plt.ylim([-0.1, 1.1])
        plt.xlim([-10, 110])
        # plt.ylim([-0.1, 1.3])
        # plt.xticks(fontsize=40)
        plt.xticks([0, 100], fontsize=20)
        # plt.yticks(fontsize=40)
        plt.yticks(np.arange(0, 1.1, 0.5), fontsize=20)
        # Remove tick lines while keeping labels
        plt.tick_params(axis='both', which='both', length=0)
        # plt.ylabel('μ', fontsize=44)
        # plt.xlabel('TCN', fontsize=40)
    # plt.legend()
    plt.show()

    # fig.savefig("Fig2DRight.pdf", bbox_inches="tight")
    return LineLabel, PCN, c, cmap, color, color1, fig, g, growthRateMatrix


@app.cell
def __(np):
    ### A function to calculate the transition matrix between neighboring subpopulations

    ### It takes 5 inputs, with last 2 being optional:
    ### plasmidNum: the number of plasmids;
    ###          If there are n plasmids, then there are n+1 maximum possible subpopulations,
    ###          each subpopulation has a different number of transposons on plasmids, from 0 to n
    ### transitionR: rate of transitions between neighboring subpopulations
    ### excisionEffect: a percentage of transitionR
    ###                 transitionR * excisionEffect is the transition rate from subpopulation n to n-1
    ###                 it's zero for incompatible plasmids

    ### It returns the transition matrix

    def calculateTransitionMatrix(plasmidNum, transitionR, excisionEffect, jumping = 1, plasmidAmp = 1):
        # For PCN dynamics effect
        ampFactor = 2

        if jumping != 1:
            # excisionEffect = 0
            transitionRForward = transitionR
            bottomValueList = transitionRForward*np.ones(plasmidNum)
            transitionRReverse = transitionR
            upValueList = transitionRReverse*np.ones(plasmidNum)
        else:
            # transposition increase the forward rate by 0.05; could adjust: 0.1 & 0.03 are presented in SF2
            transitionRForward = transitionR+0.05
            constant = 0.0;
            jumpingAsAFunctionOfTPN = np.array([constant * i for i in range(plasmidNum)])

            # jumpingAsAFunctionOfTPN = transitionRForward*np.ones(plasmidNum)

            bottomValueListBase = transitionRForward*np.ones(plasmidNum)
            bottomValueList = jumpingAsAFunctionOfTPN + bottomValueListBase

            transitionRReverse = transitionR
            upValueList = transitionRReverse*np.ones(plasmidNum)

        if plasmidAmp == 1:
            # transitionRForward = transitionRForward*1.0*2
            # transitionRReverse = transitionRReverse*2
            bottomValueList *= ampFactor
            upValueList *= ampFactor
        else:
            ampFactor=1

        # A cap for limited transposition resources
        UpBottomSum = upValueList + bottomValueList
        UpBottomSum[UpBottomSum>1] = 1
        bottomValueList[bottomValueList>(1-upValueList[2])] = 1-upValueList[2]

        # Upper part of the triD matrix, i.e., A to A-1
        # upValueList = transitionRReverse*np.ones(plasmidNum)
        upMatrix = np.diag(upValueList, 1)

        # Bottom part of the triD matrix, i.e., A to A+1
        ### This could be revised to reflect more real mechanisms
        # bottomValueList = transitionRForward*np.ones(plasmidNum)
        bottomMatrix = np.diag(bottomValueList, -1)

        # Diagonal part of the triD matrix, i.e., A's loss due to its transition to others
        diagValueList = np.zeros(plasmidNum+1) 
        diagValueList[:len(bottomValueList)] -= bottomValueList
        diagValueList[1:1+len(bottomValueList)] -= upValueList

        diagMatrix = np.diag(diagValueList)

        # Merge these three parts into one matrix for the base matrix
        finalMatrix = upMatrix + bottomMatrix + diagMatrix

        # Special consideration for the 0 -> 1 (only possible when there is transposon)
        if jumping != 1:
            finalMatrix[1, 0] = 0*finalMatrix[1, 0]
        else:
            finalMatrix[1, 0] = finalMatrix[1, 0]-ampFactor*transitionR

        finalMatrix[0, 0] = 0-finalMatrix[1, 0]

        # Special consideration for the n -> n-1 (should be very rare)
        finalMatrix[-2, -1] = excisionEffect*finalMatrix[-2, -1]
        finalMatrix[-1, -1] = 0-finalMatrix[-2, -1]

        return finalMatrix
    return calculateTransitionMatrix,


@app.cell
def __():
    ### Testing here

    # 100 plasmids, 101 subpopulations, with transition rate = 0.1, and excision effect = 0.01
    plasmidN = 100
    transR = 0.1
    excisionE = 0.0
    return excisionE, plasmidN, transR


@app.cell
def __(calculateTransitionMatrix, excisionE, plasmidN, transR):
    calculateTransitionMatrix(plasmidN, transR, excisionE, jumping = -1, plasmidAmp = -1)
    return


@app.cell
def __(calculateTransitionMatrix, excisionE, plasmidN, transR):
    calculateTransitionMatrix(plasmidN, transR, excisionE, jumping = -1, plasmidAmp = 1)
    return


@app.cell
def __(calculateTransitionMatrix, excisionE, plasmidN, transR):
    calculateTransitionMatrix(plasmidN, transR, excisionE, jumping = 1, plasmidAmp = 1)
    return


@app.cell
def __(calculateTransitionMatrix, excisionE, plasmidN, transR):
    calculateTransitionMatrix(plasmidN, transR, excisionE, jumping = 1, plasmidAmp = -1)
    return


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Here is the ODE model.
        It takes in the finalized triD matrix: growth rate + transition.
        """
    )
    return


@app.cell
def __():
    ### A function to set up the logistic eqn to be run by an ode solver

    ### It takes 4 inputs:
    ### time: the time span to run the ODE
    ### state: the initial condition, i.e. a list of population density or OD, normalized or not
    ### triDMatrix: the tridiagonal matrix, a combination of growth and transition matrix
    ### dilutionR: the periodic dilution rate, a scalar value

    ### It returns the list of differential equation values for the ODE solver to be used

    def triDLogisticGrowth(time, state, triDMatrix, dilutionR):
        # Calculate the logistic term first
        logisticTerm = 1-sum(state)

        # Now write the differential eqn
        expList = triDMatrix@state # Exponential growth term that takes both transition and growth into account
        logisticList = expList*logisticTerm # Logistic growth term
        logisticWithDiluList = logisticList - dilutionR*state # Subtract the dilution term

        return logisticWithDiluList
    return triDLogisticGrowth,


@app.cell
def __(
    calculateTransitionMatrix,
    excisionE,
    np,
    plasmidN,
    transR,
    triDLogisticGrowth,
):
    ### Testing here
    # Dilution rate
    dilutionR = 0.1

    # Initialize the y0, each of density = 0.1
    y0 = 0.1*np.ones(plasmidN+1)

    # A simple triD matrix
    # 101 subpopulations, with transition rate = 0.1, and excision effect = 0.01
    A = calculateTransitionMatrix(plasmidN, transR, excisionE)

    # Call function, t value does not matter here, as it is only used by the ODE solver later
    triDLogisticGrowth(1, y0, A, dilutionR)
    return A, dilutionR, y0


@app.cell
def __(np, pd):
    ### It is to be ran After having the simulated time course from the function sameGrowthEffectSimOnce.

    ### It takes in 2 inputs:
    ### timecourse: the simulated time course from the function sameGrowthEffectSimOnce
    ### currentGrowthRate: the list of growth rates for all subpopulations under the current atb selection
    ###                    this is from the function growthRateMatrix

    def priceResults(timecourse, currentGrowthRate, plasmidConstraint=-1):
        # The steps of the time course
        steps = timecourse.shape[1]
        # For population-level transposon copy & growth rate
        meanCopyList = []
        meanMuList = []
        # For Price THM calculation, i.e., normalized by total OD. Again, population-level
        meanCopyList_P = []
        meanMuList_P = []
        VarList = []
        CovList = []
        # Sum of OD at each time point
        sumODList = []

        # Calculate Price's THM for all time steps
        for s in np.arange(0, steps, 10):   
            # Get the ODs of all subpopulations at that time step
            yCurrent = timecourse[:,s]
            # Sum the ODs
            sumODList.append(sum(yCurrent))
            # Calculate the percentage of all subpopulations
            yCurrentPercentage = [i/sum(yCurrent) for i in yCurrent]

            # # Get the logistic term for the calculation of effective growth rates
            logisticTerm = 1-sum(yCurrent)

            # # To record mean transposon copy number and mean effective growth rate
            # meanCopy = 0
            # meanMu = 0
            meanCopy_P = 0
            meanMu_P = 0

            # # Calculate mean transposon copy number and mean eff growth rate
            for c in range(len(yCurrent)):
                if plasmidConstraint == -1:
                    currCopy = c
                else:
                    currCopy = c+len(yCurrent)-1
                # Calculate normalized mean copy number using the percentage
                currCopyMean_P = yCurrentPercentage[c]*currCopy 
                meanCopy_P = meanCopy_P + currCopyMean_P

                # Calculate normalized mean effective growth rate using the percentage
                currMuMean_P = yCurrentPercentage[c]*currentGrowthRate[c]*logisticTerm
                meanMu_P = meanMu_P + currMuMean_P

            # meanCopyList.append(meanCopy)
            # meanMuList.append(meanMu)
            meanCopyList_P.append(meanCopy_P)
            meanMuList_P.append(meanMu_P)


        ### Return as a Dataframe for pythonic calculation
        timeList = list(np.arange(0, len(meanCopyList_P)))

        df = pd.DataFrame(list(zip(timeList, sumODList, meanCopyList_P, meanMuList_P)), \
                          columns=['Time', 'OD', 'TCN', 'meanMu'])
        return df
    return priceResults,


@app.cell
def __(np):
    ### It takes one input:
    ### meanCopyList: The list of mean transposon copies over the time course
    ### It returns the response speed.

    def responseSpeed(meanCopyList, ODList):     
        # Calculate all transoposon copy changing rate over time   
        meanCopyFinal = meanCopyList[-1]

        # if meanCopyFinal < 0.01: # effectively 0, for division purpose, give it a small value
        #     meanCopyFinal = 0.01

        if meanCopyList[0] != 0:
            meanCopyInitial = meanCopyList[0]
        else:
            meanCopyInitial = meanCopyList[1]

        numberGeneration = np.log2(ODList[-1]/ODList[0])

        if ODList[-1] > 0.1:
            rate = np.log(meanCopyFinal/meanCopyInitial)/numberGeneration
        else:
            rate = 0


        return rate
    return responseSpeed,


@app.cell
def __(
    calculateTransitionMatrix,
    np,
    priceResults,
    responseSpeed,
    solve_ivp,
    stats,
    t_step,
    triDLogisticGrowth,
):
    ### A function to simulate the population under one atb concentration.

    ### It takes 7-11 inputs (the last 4 are optional):
    ###
    ### The mandatory inputs:
    ### plasmidNum: the number of plasmids;
    ###          If there are n plasmids, then there are n+1 maximum possible subpopulations,
    ###          each subpopulation has a different number of transposons on plasmids, from 0 to n
    ### transitionR: rate of transitions between neighboring subpopulations
    ### excisionEffect: the percentage of diffusionR
    ###                 diffusionR * excisionEffect is the transition rate from subpopulation n to n-1
    ###                 & the transition rate from subpopulation 0 to 1
    ### dilutionR: the periodic dilution rate, a scalar value
    ### time: the time length to simulate the ODE system for
    ### growthRate: the full list of growth rates of a population with max possible plasmids,
    ###             under a given selection

    ####### + the plasmidNum here to a range, take the range of the baslineGrowthRate - > ratio is fine
    def SimOnce(plasmidNum, transitionR, excisionEffect, dilutionR, time, growthRate, baselineGrowthRate, \
                jumping = 1, plasmidAmp = 1, y0=-1, Up=1):
        # If no valid y0 list at t0, initialize it
        popSize = int(plasmidNum+1)

        if np.size(y0) == 1:
            if Up == -1: # no selection
                # Normal distribution
                x = np.arange(0, popSize)
                pdf = stats.norm.pdf(x, 0.80*popSize, 0.085*popSize) # both values can vary
                # pdf = stats.norm.pdf(x, 0.90*popSize, 0.30*popSize) # for exhaustive comparison of response speed over PCN; Fig 5
                y0 = pdf*0.000000001

            else: # with selection
                x = np.arange(0, popSize)
                pdf = stats.norm.pdf(x, 0.20*popSize, 0.085*popSize) # both values can vary
                # pdf = stats.norm.pdf(x, 10, 0.30*popSize) # for exhaustive comparison of response speed over PCN; Fig 5
                y0 = pdf*0.000000001


        # Calculate the transition matrix first
        transitionMatrix = calculateTransitionMatrix(plasmidNum, transitionR, excisionEffect, jumping, plasmidAmp)

        # Combine the transition with growth rate for the final triD matrix
        finalMatrix = transitionMatrix + np.diag(growthRate)

        # Set up to run the ODE
        # Run time
        t_span = (0.0, time) 

        parameters = (finalMatrix, dilutionR)       
        # # Run ODE solver
        result_solve_ivp = solve_ivp(triDLogisticGrowth, t_span, y0,
                                     method='LSODA', args=parameters, \
                                     max_step = t_step, atol = 1, rtol = 1)

        ## Figures for trouble shooting

        # Plot the time course
        # fig = plt.figure(figsize=(5, 3))
        # for species in range(result_solve_ivp.y.shape[0]):
        #     speciesLabel = 'copy number = ' + str(species)
        #     plt.plot(result_solve_ivp.t, result_solve_ivp.y[species, :], label=speciesLabel)
        # # Only print the species number and color legend if there are < 15 species
        # if species < 15:
        #     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # plt.ylabel('OD', fontsize=20)
        # plt.yticks(fontsize=20)
        # plt.xlabel('Time') # , fontsize=20
        # plt.xticks(fontsize=20)
        # plt.title('Time Course', fontsize = 25)
        # plt.show()

        # fig = plt.figure(figsize=(5, 3))
        # for species in range(result_solve_ivp.y.shape[0]):
        #     speciesLabel = 'copy number = ' + str(species)
        #     plt.semilogy(result_solve_ivp.t, result_solve_ivp.y[species, :], label=speciesLabel)
        # # Only print the species number and color legend if there are < 15 species
        # if species < 15:
        #     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # plt.ylabel('OD', fontsize=20)
        # plt.yticks(fontsize=20)
        # plt.xlabel('Time') # , fontsize=20
        # plt.xticks(fontsize=20)
        # plt.title('Time Course', fontsize = 25)
        # plt.show()

        # First time point's ys
        # fig = plt.figure(figsize=(5, 3)) 
        # plt.bar(range(len(result_solve_ivp.y[:, 0])), result_solve_ivp.y[:, 0])
        # plt.ylabel('OD', fontsize=20)
        # plt.yticks(fontsize=20)
        # plt.xlabel('Copy Number') #, fontsize=20
        # plt.xticks(fontsize=20)
        # plt.title('Initial Time Point', fontsize = 25)    
        # plt.show()

        # Last time point's ys
        # fig = plt.figure(figsize=(5, 3)) 
        # plt.bar(range(len(result_solve_ivp.y[:, -1])), result_solve_ivp.y[:, -1])
        # plt.ylabel('OD', fontsize=20)
        # plt.yticks(fontsize=20)
        # plt.xlabel('Copy Number') #, fontsize=20
        # plt.xticks(fontsize=20)
        # plt.title('Final Time Point', fontsize = 25)    
        # plt.show()

        # The full time courses for all subpopulations
        simulatedY = result_solve_ivp.y
        simulatedTime = result_solve_ivp.t

        SimOnce_DF = \
        priceResults(simulatedY, growthRate)

        # Store the response speed
        Response_Speed = responseSpeed(SimOnce_DF['TCN'].values, SimOnce_DF['OD'].values)
        if Up == -1:
            print(Response_Speed)
            Response_Speed = -1*(Response_Speed)
            print(Response_Speed)


        SimOnce_DF['response speed'] = Response_Speed
        SimOnce_DF['transition'] = transitionR
        SimOnce_DF['jumping'] = jumping
        SimOnce_DF['plasmidAmp'] = plasmidAmp

        yEnd = simulatedY[:,-1]

        return SimOnce_DF, yEnd #, result_solve_ivp
    return SimOnce,


@app.cell
def __():
    return


@app.cell
def __(SimOnce, growthRateFluctuating, plt):
    #Test here
    muBaseline = 1

    # Caculate the growth rate for up now
    growthMatrix = growthRateFluctuating(100, muBaseline)
    currentGrowthRateUp = growthMatrix.T[-1]
    currentGrowthRateDown = growthMatrix.T[0]

    fig02 = plt.figure(figsize=(2,2))
    plt.plot(currentGrowthRateUp, linewidth = 3)
    plt.title('growth rates under selection', fontsize = 20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()
    # Run once given the current growth rate set above
    transitionR_test = 0.1
    excisionE_test = 1.0
    dilutionR_test = 0.2
    time = 2000

    # currentGrowthRate_test, simulatedY_test, simulatedTime_test = \
    SimOnce_DF_test, yEnd_test = \
    SimOnce(100, transitionR_test, excisionE_test, dilutionR_test, time, currentGrowthRateUp, muBaseline, 1)
    return (
        SimOnce_DF_test,
        currentGrowthRateDown,
        currentGrowthRateUp,
        dilutionR_test,
        excisionE_test,
        fig02,
        growthMatrix,
        muBaseline,
        time,
        transitionR_test,
        yEnd_test,
    )


@app.cell
def __(SimOnce, growthRateFluctuating, np, pd):
    ### Add in the growth rate to be normalized

    def SimOnceMultiplePCN(muBaselineList, plasmidNumList, transitionR, excisionE, dilutionR, time, \
                           jumping = 1, plasmidAmp = 1, y0 = -1, Up = 1):   
        finalDF = pd.DataFrame()
        yEndFinalDF = pd.DataFrame()

        # Now for the loop for each PCN
        for j in range(len(plasmidNumList)):
            p = plasmidNumList[j]
            muBaseline = muBaselineList[j]
            if len(transitionR) == 1:
                transR = transitionR[0]
            else:
                transR = transitionR[j]

            growthMatrix = growthRateFluctuating(p, muBaseline)

            if Up == 1:
                currentGrowthRate = growthMatrix.T[-1]
            else:
                # Calculate the growth rate for down
                currentGrowthRate = growthMatrix.T[0]   

            # print(y0)

            # If y0 is predetermined
            if isinstance(y0, pd.DataFrame):
                y0Use = y0[y0['PCN']==p]['lastY']      
            elif y0 == -1:
                y0Use = y0


            Once_DF, yEnd = SimOnce(p, transR, excisionE, dilutionR, time, currentGrowthRate, muBaseline, \
                              jumping, plasmidAmp, y0Use, Up) # np.max(currentGrowthRate)-np.min(currentGrowthRate)
            Once_DF['PCN'] = p    
            finalDF = pd.concat([finalDF, Once_DF], ignore_index=True)

            yEndDF = pd.DataFrame({
                'lastY': yEnd,
                'PCN': p
            })
            yEndFinalDF = pd.concat([yEndFinalDF, yEndDF], ignore_index=True)

        if Up == 1:
            finalDF['Selection'] = 1
        else:
            finalDF['Selection'] = -1

        finalDF['log10(TCN)'] = np.log10(finalDF['TCN'])

        # print(finalDF)
        return finalDF, yEndFinalDF
    return SimOnceMultiplePCN,


@app.cell
def __(mo):
    mo.md(r"# Fig. 2")
    return


@app.cell
def __(muBaselineList):
    # Initialize all parameters
    muBaselineList_one = [muBaselineList[-1]]
    excisionE_ = 0.1
    dilutionR_ = 0.2
    time_S = 2000
    time_NoS = 3000
    y0_ = -1
    Up = -1
    plasmidN_ = 100
    return (
        Up,
        dilutionR_,
        excisionE_,
        muBaselineList_one,
        plasmidN_,
        time_NoS,
        time_S,
        y0_,
    )


@app.cell
def __(
    SimOnceMultiplePCN,
    Up,
    excisionE_,
    muBaselineList_one,
    plasmidN_,
    time_S,
    y0_,
):
    Once_LessTransition_Selection_DF, Once_LessTransition_Selection_yFinal = \
        SimOnceMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], excisionE_, 0.15, time_S, -1, -1, y0_, (0-Up))

    Once_WithTransition_Selection_DF, Once_WithTransition_Selection_yFinal = \
        SimOnceMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], excisionE_, 0.15, time_S, 1, -1, y0_, (0-Up))

    Once_WithTransition_Selection_LargePlasmidRange_DF, Once_WithTransition_Selection_LargePlasmidRange_yFinal = \
        SimOnceMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], excisionE_, 0.15, time_S, 1, 1, y0_, (0-Up))
    return (
        Once_LessTransition_Selection_DF,
        Once_LessTransition_Selection_yFinal,
        Once_WithTransition_Selection_DF,
        Once_WithTransition_Selection_LargePlasmidRange_DF,
        Once_WithTransition_Selection_LargePlasmidRange_yFinal,
        Once_WithTransition_Selection_yFinal,
    )


@app.cell
def __(
    SimOnceMultiplePCN,
    Up,
    dilutionR_,
    excisionE_,
    muBaselineList_one,
    plasmidN_,
    time_NoS,
    y0_,
):
    Once_LessTransition_NoSelection_DF, Once_LessTransition_NoSelection_yFinal = \
        SimOnceMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], excisionE_, dilutionR_, time_NoS, -1, -1, y0_, Up)

    Once_WithTransition_NoSelection_DF, Once_WithTransition_NoSelection_yFinal = \
        SimOnceMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], excisionE_, dilutionR_, time_NoS, 1, -1, y0_, Up)

    Once_WithTransition_NoSelection_LargePlasmidRange_DF, Once_WithTransition_NoSelection_LargePlasmidRange_yFinal = \
        SimOnceMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], excisionE_, dilutionR_, time_NoS, 1, 1, y0_, Up)
    return (
        Once_LessTransition_NoSelection_DF,
        Once_LessTransition_NoSelection_yFinal,
        Once_WithTransition_NoSelection_DF,
        Once_WithTransition_NoSelection_LargePlasmidRange_DF,
        Once_WithTransition_NoSelection_LargePlasmidRange_yFinal,
        Once_WithTransition_NoSelection_yFinal,
    )


@app.cell
def __(
    Once_LessTransition_NoSelection_DF,
    Once_LessTransition_Selection_DF,
    Once_WithTransition_NoSelection_DF,
    Once_WithTransition_NoSelection_LargePlasmidRange_DF,
    Once_WithTransition_Selection_DF,
    Once_WithTransition_Selection_LargePlasmidRange_DF,
):
    Once_LessTransition_NoSelection_DF['κb;κf'] = '0.1; 0.1'
    Once_WithTransition_NoSelection_DF['κb;κf'] = '0.1; 0.1+'
    Once_WithTransition_NoSelection_LargePlasmidRange_DF['κb;κf'] = '0.2; 0.2+'

    Once_LessTransition_Selection_DF['κb;κf'] = '0.1; 0.1'
    Once_WithTransition_Selection_DF['κb;κf'] = '0.1; 0.1+'
    Once_WithTransition_Selection_LargePlasmidRange_DF['κb;κf'] = '0.2; 0.2+'
    return


@app.cell
def __(
    Once_LessTransition_NoSelection_DF,
    Once_LessTransition_Selection_DF,
    Once_WithTransition_NoSelection_DF,
    Once_WithTransition_NoSelection_LargePlasmidRange_DF,
    Once_WithTransition_Selection_DF,
    Once_WithTransition_Selection_LargePlasmidRange_DF,
    pd,
):
    results_selection = pd.concat([Once_LessTransition_Selection_DF, Once_WithTransition_Selection_DF, Once_WithTransition_Selection_LargePlasmidRange_DF])

    results_noS = pd.concat([Once_LessTransition_NoSelection_DF, Once_WithTransition_NoSelection_DF, Once_WithTransition_NoSelection_LargePlasmidRange_DF])
    return results_noS, results_selection


@app.cell
def __(alt, plasmidN_, results_noS):
    # y axis is log(TCN)
    F2ELeft = alt.Chart(results_noS).mark_line(
        point=alt.OverlayMarkDef(size=10, opacity=0.7)).encode(
        x=alt.X('Time', scale=alt.Scale(domain=[0, 1800]), axis=alt.Axis(labelFontSize=30, labelAngle=0, ticks=False, grid=False, title= None)),
        y=alt.Y('TCN', scale=alt.Scale(type="log", domain=[0.01, plasmidN_]), axis=alt.Axis(labelFontSize=30, ticks=False, values=[0.01, 1, plasmidN_], grid=False, title= None)), #, scale=alt.Scale(type="log", domain=[0.01, 100])
        color=alt.Color('κb;κf:O', legend=None).scale(scheme="purples", reverse=False), #greens, goldgreen, warmgreys
    ).properties(
        width=150,
        height=150
    ).interactive()

    F2ELeft
    return F2ELeft,


@app.cell
def __(F2ELeft):
    F2ELeft.save('SF3B_muDown.pdf')
    return


@app.cell
def __(alt, results_noS):
    # y axis is TCN
    alt.Chart(results_noS).mark_line(
        point=alt.OverlayMarkDef(size=10, opacity=0.7)).encode(
        x=alt.X('Time', scale=alt.Scale(domain=[0, 1800]), axis=alt.Axis(labelFontSize=30, labelAngle=0, ticks=False, grid=False, title= None)),
        y=alt.Y('TCN', scale=alt.Scale(domain=[0, 100]), axis=alt.Axis(labelFontSize=30, ticks=False, tickCount=2, grid=False, title= None)), #, scale=alt.Scale(type="log", domain=[0.01, 100])
        color=alt.Color('κb;κf:O', legend=None).scale(scheme="purples", reverse=False), #greens, goldgreen, warmgreys
    ).properties(
        width=150,
        height=150
    ).interactive()
    return


@app.cell
def __():
    return


@app.cell
def __(alt, plasmidN_, results_selection):
    # y axis is log(TCN)
    F2ERight = alt.Chart(results_selection[(results_selection['κb;κf']=='0.1; 0.1') | (results_selection['κb;κf']=='0.1; 0.1+') | (results_selection['κb;κf']=='0.2; 0.2+')]).mark_line(
        point=alt.OverlayMarkDef(size=10, opacity=0.7)).encode(
        x=alt.X('Time', scale=alt.Scale(domain=[0, 500]), axis=alt.Axis(labelFontSize=30, labelAngle=0, ticks=False, tickCount=1, grid=False, title= None)),
        y=alt.Y('TCN', scale=alt.Scale(domain=[10, plasmidN_], type="log"), axis=alt.Axis(labelFontSize=30, ticks=False, values=[10, 20, plasmidN_], grid=False, title= None)),
        color=alt.Color('κb;κf:O', legend=None, scale=alt.Scale(domain=['0.1; 0.1', '0.1; 0.1+', '0.2; 0.2+'], range=['#C3C5DE', '#A2A0C8', '#644F9C']))
    ).properties(
        width=150,
        height=150
    ).interactive()

    F2ERight
    return F2ERight,


@app.cell
def __(F2ERight):
    F2ERight.save('SF3B_muUp.pdf')
    return


@app.cell
def __(alt, results_selection):
    # y axis is TCN
    alt.Chart(results_selection).mark_line(
        point=alt.OverlayMarkDef(size=10, opacity=0.7)).encode(
        x=alt.X('Time', scale=alt.Scale(domain=[0, 500]), axis=alt.Axis(labelFontSize=30, labelAngle=0, ticks=False, tickCount=1, grid=False, title= None)),
        y=alt.Y('TCN', scale=alt.Scale(domain=[0, 100]), axis=alt.Axis(labelFontSize=30, ticks=False, tickCount=2, grid=False, title= None)), #, scale=alt.Scale(domain=[2.6, 5])
        color=alt.Color('κb;κf:O', legend=None).scale(scheme="purples", reverse=False), #greens, goldgreen, warmgreys
    ).properties(
        width=150,
        height=150
    ).interactive()
    return


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"# Fig.5")
    return


@app.cell
def __(alt):
    alt.data_transformers.disable_max_rows()
    return


@app.cell
def __(SimOnceMultiplePCN, Up, muBaselineList, popSizeList, y0_):
    # Exhaustive, no selection

    Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_DF_02, Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_yFinal_02 = \
        SimOnceMultiplePCN(muBaselineList, popSizeList, [0.2], 0.0, 0.1, 5000, -1, -1, y0_, Up)

    Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_DF_02, Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_yFinal_02 = \
        SimOnceMultiplePCN(muBaselineList, popSizeList, [0.2], 0.0, 0.1, 5000, -1, 1, y0_, Up)

    Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_DF_02, Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_yFinal_02 = \
        SimOnceMultiplePCN(muBaselineList, popSizeList, [0.2], 0.0, 0.1, 5000, 1, -1, y0_, Up)

    Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_DF_02, Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_yFinal_02 = \
        SimOnceMultiplePCN(muBaselineList, popSizeList, [0.2], 0.0, 0.1, 5000, 1, 1, y0_, Up)


    Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_ResponseSpeed_02 = Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_DF_02[['response speed', 'PCN']].drop_duplicates()

    Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_ResponseSpeed_02 = Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_DF_02[['response speed', 'PCN']].drop_duplicates()

    Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_ResponseSpeed_02 = Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_DF_02[['response speed', 'PCN']].drop_duplicates()

    Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_ResponseSpeed_02 = Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_DF_02[['response speed', 'PCN']].drop_duplicates()
    return (
        Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_DF_02,
        Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_ResponseSpeed_02,
        Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_yFinal_02,
        Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_DF_02,
        Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_ResponseSpeed_02,
        Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_yFinal_02,
        Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_DF_02,
        Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_ResponseSpeed_02,
        Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_yFinal_02,
        Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_DF_02,
        Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_ResponseSpeed_02,
        Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_yFinal_02,
    )


@app.cell
def __(
    Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_ResponseSpeed_02,
    alt,
    figureHeight,
    figureWidth,
):
    # without selection

    F5AB = speed_LessT_NoS_LessPlasmidRange_Exhaustive_02 = alt.Chart(Once_LessTransition_NoSelection_LessPlasmidRange_Exhaustive_ResponseSpeed_02).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('PCN', scale=alt.Scale(domain=[10, 100]), axis=alt.Axis(labelFontSize=30, values=[10, 100], labelAngle=0, ticks=False, grid=False, title=None)),
        y=alt.Y('response speed', scale=alt.Scale(domain=[-0.001, 0.4]),  axis=alt.Axis(labelFontSize=30, tickCount=2, ticks=False, grid=False, title=None)), #, scale=alt.Scale(domain=[0, 0.255]), 
        color=alt.Color('PCN:O', legend=None,).scale(scheme="plasma",reverse=True),
    ).properties(
        width=figureWidth,
        height=figureHeight
    ).interactive()


    speed_LessT_NoS_LessPlasmidRange_Exhaustive_02.configure_axis(
        labelFontSize=20,
        titleFontSize=20
    ).configure_header(
        titleFontSize=20,
        labelFontSize=20 
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )

    F5AB
    return F5AB, speed_LessT_NoS_LessPlasmidRange_Exhaustive_02


@app.cell
def __():
    # F5AB.save('SF9A_mean50_std20_noT.pdf')
    return


@app.cell
def __(
    Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_ResponseSpeed_02,
    alt,
    figureHeight,
    figureWidth,
):
    # without selection

    speed_LessT_NoS_LargePlasmidRange_Exhaustive_02 = alt.Chart(Once_LessTransition_NoSelection_LargePlasmidRange_Exhaustive_ResponseSpeed_02).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('PCN', scale=alt.Scale(domain=[10, 100]), axis=alt.Axis(labelFontSize=30, values=[10, 100], labelAngle=0, ticks=False, grid=False, title=None)),
        y=alt.Y('response speed', scale=alt.Scale(domain=[-0.001, 0.45]),  axis=alt.Axis(labelFontSize=30, tickCount=2, ticks=False, grid=False, title=None)), #, scale=alt.Scale(domain=[0, 0.255]), 
        color=alt.Color('PCN:O', legend=None,).scale(scheme="plasma",reverse=True),
    ).properties(
        width=figureWidth,
        height=figureHeight
    ).interactive()


    speed_LessT_NoS_LargePlasmidRange_Exhaustive_02.configure_axis(
        labelFontSize=20,
        titleFontSize=20
    ).configure_header(
        titleFontSize=20,
        labelFontSize=20 
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )
    return speed_LessT_NoS_LargePlasmidRange_Exhaustive_02,


@app.cell
def __(
    Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_ResponseSpeed_02,
    alt,
    figureHeight,
    figureWidth,
):
    # without selection

    F5AA = speed_WithT_NoS_LessPlasmidRange_Exhaustive_02 = alt.Chart(Once_WithTransition_NoSelection_LessPlasmidRange_Exhaustive_ResponseSpeed_02).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('PCN', scale=alt.Scale(domain=[10, 100]), axis=alt.Axis(labelFontSize=30, values=[10, 100], labelAngle=0, ticks=False, grid=False, title=None)),
        y=alt.Y('response speed', scale=alt.Scale(domain=[-0.001, 0.2]),  axis=alt.Axis(labelFontSize=30, tickCount=2, ticks=False, grid=False, title=None)), #, scale=alt.Scale(domain=[0, 0.255]), 
        color=alt.Color('PCN:O', legend=None,).scale(scheme="plasma",reverse=True),
    ).properties(
        width=figureWidth,
        height=figureHeight
    ).interactive()


    speed_WithT_NoS_LessPlasmidRange_Exhaustive_02.configure_axis(
        labelFontSize=20,
        titleFontSize=20
    ).configure_header(
        titleFontSize=20,
        labelFontSize=20 
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )

    F5AA
    return F5AA, speed_WithT_NoS_LessPlasmidRange_Exhaustive_02


@app.cell
def __():
    # F5AA.save('SF9A_mean50_std20_T.pdf')
    return


@app.cell
def __(
    Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_ResponseSpeed_02,
    alt,
    figureHeight,
    figureWidth,
):
    # without selection

    speed_WithT_NoS_LargePlasmidRange_Exhaustive_02 = alt.Chart(Once_WithTransition_NoSelection_LargePlasmidRange_Exhaustive_ResponseSpeed_02).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('PCN', scale=alt.Scale(domain=[10, 100]), axis=alt.Axis(labelFontSize=30, values=[10, 100], labelAngle=0, ticks=False, grid=False, title=None)),
        y=alt.Y('response speed', scale=alt.Scale(domain=[-0.001, 0.2]),  axis=alt.Axis(labelFontSize=30, tickCount=2, ticks=False, grid=False, title=None)), #, scale=alt.Scale(domain=[0, 0.255]), 
        color=alt.Color('PCN:O', legend=None,).scale(scheme="plasma",reverse=True),
    ).properties(
        width=figureWidth,
        height=figureHeight
    ).interactive()


    speed_WithT_NoS_LargePlasmidRange_Exhaustive_02.configure_axis(
        labelFontSize=20,
        titleFontSize=20
    ).configure_header(
        titleFontSize=20,
        labelFontSize=20 
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )
    return speed_WithT_NoS_LargePlasmidRange_Exhaustive_02,


@app.cell
def __():
    return


@app.cell
def __(SimOnceMultiplePCN, Up, muBaselineList, popSizeList, y0_):
    # Exhaustive, selection

    Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_DF_02, Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_yFinal_02 = \
        SimOnceMultiplePCN(muBaselineList, popSizeList, [0.4], 0.1, 0.1, 5000, -1, -1, y0_, (0-Up))

    Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_DF_02, Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_yFinal_02 = \
        SimOnceMultiplePCN(muBaselineList, popSizeList, [0.4], 0.1, 0.1, 5000, -1, 1, y0_, (0-Up))

    Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_DF_02, Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_yFinal_02 = \
        SimOnceMultiplePCN(muBaselineList, popSizeList, [0.4], 0.1, 0.1, 5000, 1, -1, y0_, (0-Up))

    Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_DF_02, Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_yFinal_02 = \
        SimOnceMultiplePCN(muBaselineList, popSizeList, [0.4], 0.1, 0.1, 5000, 1, 1, y0_, (0-Up))

    Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_ResponseSpeed_02 = Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_DF_02[['response speed', 'PCN']].drop_duplicates()

    Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_ResponseSpeed_02 = Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_DF_02[['response speed', 'PCN']].drop_duplicates()

    Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_ResponseSpeed_02 = Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_DF_02[['response speed', 'PCN']].drop_duplicates()

    Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_ResponseSpeed_02 = Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_DF_02[['response speed', 'PCN']].drop_duplicates()
    return (
        Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_DF_02,
        Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_ResponseSpeed_02,
        Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_yFinal_02,
        Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_DF_02,
        Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_ResponseSpeed_02,
        Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_yFinal_02,
        Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_DF_02,
        Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_ResponseSpeed_02,
        Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_yFinal_02,
        Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_DF_02,
        Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_ResponseSpeed_02,
        Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_yFinal_02,
    )


@app.cell
def __(
    Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_ResponseSpeed_02,
    alt,
    figureHeight,
    figureWidth,
):
    # under selection

    F5AD = speed_LessTran_S_LessPlasmidRange_Exhaustive_02 = alt.Chart(Once_LessTransition_Selection_LessPlasmidRange_Exhaustive_ResponseSpeed_02).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('PCN', axis=alt.Axis(values=[10, 100], labelFontSize=30, labelAngle=0, ticks=False, grid=False, title = None)),
        y=alt.Y('response speed', scale=alt.Scale(domain=[-0.001, 0.045]), axis=alt.Axis(labelFontSize=30, ticks=False, tickCount=2, grid=False, title = None)), #, scale=alt.Scale(domain=[2.6, 5]), , scale=alt.Scale(domain=[0, 1.6]), domain=[0, 0.22]
        color=alt.Color('PCN:O', legend=None).scale(scheme="plasma",reverse=True),
    ).properties(
        width=figureWidth,
        height=figureHeight
    ).interactive()


    speed_LessTran_S_LessPlasmidRange_Exhaustive_02.configure_axis(
        labelFontSize=20,
        titleFontSize=20
    ).configure_header(
        titleFontSize=20,
        labelFontSize=20 
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )

    F5AD
    return F5AD, speed_LessTran_S_LessPlasmidRange_Exhaustive_02


@app.cell
def __(F5AD):
    # F5AD.save('SF9B_mean5_std20_noT.pdf')
    F5AD.save('SF10C_mean10_std30_noT_withK.pdf')
    return


@app.cell
def __(
    Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_ResponseSpeed_02,
    alt,
    figureHeight,
    figureWidth,
):
    SF10_noT_noK = speed_LessT_S_LargePlasmidRange_Exhaustive_02 = alt.Chart(Once_LessTransition_Selection_LargePlasmidRange_Exhaustive_ResponseSpeed_02).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('PCN', axis=alt.Axis(values=[10, 100], labelFontSize=30, labelAngle=0, ticks=False, grid=False, title = None)),
        y=alt.Y('response speed', scale=alt.Scale(domain=[-0.001, 0.045]), axis=alt.Axis(labelFontSize=30, ticks=False, tickCount=2, grid=False, title = None)), #, scale=alt.Scale(domain=[2.6, 5]), , scale=alt.Scale(domain=[0, 1.6]), domain=[0, 0.22]
        color=alt.Color('PCN:O', legend=None).scale(scheme="plasma",reverse=True),
    ).properties(
        width=figureWidth,
        height=figureHeight
    ).interactive()


    speed_LessT_S_LargePlasmidRange_Exhaustive_02.configure_axis(
        labelFontSize=20,
        titleFontSize=20
    ).configure_header(
        titleFontSize=20,
        labelFontSize=20 
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )

    SF10_noT_noK
    return SF10_noT_noK, speed_LessT_S_LargePlasmidRange_Exhaustive_02


@app.cell
def __():
    # SF10_noT_noK.save('SF10C_mean10_std30_noT_noK.pdf')
    return


@app.cell
def __(
    Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_ResponseSpeed_02,
    alt,
    figureHeight,
    figureWidth,
):
    # under selection

    F5AC = speed_WithT_S_LessPlasmidRange_Exhaustive_02 = alt.Chart(Once_WithTransition_Selection_LessPlasmidRange_Exhaustive_ResponseSpeed_02).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('PCN', axis=alt.Axis(values=[10, 100], labelFontSize=30, labelAngle=0, ticks=False, grid=False, title = None)),
        y=alt.Y('response speed', scale=alt.Scale(domain=[-0.001, 0.045]), axis=alt.Axis(labelFontSize=30, tickCount=2, ticks=False, grid=False, title = None)), #, scale=alt.Scale(domain=[2.6, 5]), , scale=alt.Scale(domain=[0, 1.6]), domain=[0, 0.22]
        color=alt.Color('PCN:O', legend=None).scale(scheme="plasma",reverse=True),
    ).properties(
        width=figureWidth,
        height=figureHeight
    ).interactive()


    speed_WithT_S_LessPlasmidRange_Exhaustive_02.configure_axis(
        labelFontSize=20,
        titleFontSize=20
    ).configure_header(
        titleFontSize=20,
        labelFontSize=20 
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )

    F5AC
    return F5AC, speed_WithT_S_LessPlasmidRange_Exhaustive_02


@app.cell
def __():
    # F5AC.save('SF9B_mean5_std20_T.pdf')
    # F5AC.save('SF10C_mean10_std30_withT_withK.pdf')
    return


@app.cell
def __(
    Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_ResponseSpeed_02,
    alt,
    figureHeight,
    figureWidth,
):
    # under selection
    SF10_T_noK = speed_WithT_S_LargePlasmidRange_Exhaustive_02 = alt.Chart(Once_WithTransition_Selection_LargePlasmidRange_Exhaustive_ResponseSpeed_02).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('PCN', axis=alt.Axis(values=[10, 100], labelFontSize=30, labelAngle=0, ticks=False, grid=False, title = None)),
        y=alt.Y('response speed', scale=alt.Scale(domain=[-0.001, 0.045]), axis=alt.Axis(labelFontSize=30, tickCount=2, ticks=False, grid=False, title = None)), #, scale=alt.Scale(domain=[2.6, 5]), , scale=alt.Scale(domain=[0, 1.6]), domain=[0, 0.22]
        color=alt.Color('PCN:O', legend=None).scale(scheme="plasma",reverse=True),
    ).properties(
        width=figureWidth,
        height=figureHeight
    ).interactive()


    speed_WithT_S_LargePlasmidRange_Exhaustive_02.configure_axis(
        labelFontSize=20,
        titleFontSize=20
    ).configure_header(
        titleFontSize=20,
        labelFontSize=20 
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )

    SF10_T_noK
    return SF10_T_noK, speed_WithT_S_LargePlasmidRange_Exhaustive_02


@app.cell
def __():
    # SF10_T_noK.save('SF10C_mean10_std30_withT_noK.pdf')
    return


@app.cell
def __(mo):
    mo.md(r"# Fluctuating env simulation")
    return


@app.cell
def __(SimOnceMultiplePCN, np, pd):
    ### Add in the growth rate to be normalized

    def SimCycleMultiplePCN(muBaselineList, plasmidNumList, transitionR, excisionE, dilutionR, timePair, \
                           jumping = 1, plasmidAmp = 1, y0 = -1, Up = -1, cycle = 3):   
        finalCycleDF = pd.DataFrame()

        # Now for each cycle
        for c in range(cycle):
            oneCycleDF = pd.DataFrame()
            # Going down once
            Once_DF_Down, YEnd_Down = SimOnceMultiplePCN(muBaselineList, plasmidNumList, transitionR, excisionE, dilutionR, timePair[0], jumping, plasmidAmp, y0, Up)
            Once_DF_Down = Once_DF_Down[['Time', 'OD', 'TCN', 'PCN', 'transition', 'jumping', 'plasmidAmp', 'Selection', 'log10(TCN)']]
            # print(list(Once_DF_Down['TCN']))
            print('TCN: ', list(Once_DF_Down['TCN'])[-1])
            print('OD: ', list(Once_DF_Down['OD'])[-1])

            if list(Once_DF_Down['TCN'])[-1] < 1.1:
                oneCycleDF = Once_DF_Down
                finalCycleDF = pd.concat([finalCycleDF, oneCycleDF], axis=0, ignore_index=True)
                break
            if np.isnan(list(Once_DF_Down['TCN'])[-1]) == True:
                oneCycleDF = Once_DF_Down
                finalCycleDF = pd.concat([finalCycleDF, oneCycleDF], axis=0, ignore_index=True)
                break

            ## Get into the YEnd_Down YEnd column and dilute
            YEnd_Down['lastY'] = YEnd_Down['lastY'] / 1000000000
            ### Give next round the DF format & call using the PCN

            # Going up once
            Once_DF_Up, YEnd_Up = SimOnceMultiplePCN(muBaselineList, plasmidNumList, transitionR, excisionE, dilutionR, timePair[1], jumping, plasmidAmp, YEnd_Down, -1*Up)
            Once_DF_Up = Once_DF_Up[['Time', 'OD', 'TCN', 'PCN', 'transition', 'jumping', 'plasmidAmp', 'Selection', 'log10(TCN)']]
            # Concate the down and up dataframes
            oneCycleDF = pd.concat([Once_DF_Down, Once_DF_Up], axis=0, ignore_index=True)

            YEnd_Up['lastY'] = YEnd_Up['lastY'] / 1000000000
            y0 = YEnd_Up

            # Concate for all cycles
            finalCycleDF = pd.concat([finalCycleDF, oneCycleDF], axis=0, ignore_index=True)
            # print('once end')

        finalCycleDF = finalCycleDF.drop('Time', axis=1)
        finalCycleDF = finalCycleDF.reset_index(names='Time')

        # print(finalDF)
        return finalCycleDF
    return SimCycleMultiplePCN,


@app.cell
def __():
    DownT = 1400
    UpT = 1100
    return DownT, UpT


@app.cell
def __(
    DownT,
    SimCycleMultiplePCN,
    Up,
    UpT,
    dilutionR_,
    muBaselineList_one,
    plasmidN_,
    y0_,
):
    Oscilation_DF_Baseline = \
        SimCycleMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], 1, dilutionR_, [DownT, UpT], -1, -1, y0_, Up, 3)

    Oscilation_DF_Transposition = \
        SimCycleMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], 1, dilutionR_, [DownT, UpT], 1, -1, y0_, Up, 3)

    Oscilation_DF_TPV = \
        SimCycleMultiplePCN(muBaselineList_one, [plasmidN_], [0.1], 1, dilutionR_, [DownT, UpT], 1, 1, y0_, Up, 3)
    return (
        Oscilation_DF_Baseline,
        Oscilation_DF_TPV,
        Oscilation_DF_Transposition,
    )


@app.cell
def __(
    Oscilation_DF_Baseline,
    Oscilation_DF_TPV,
    Oscilation_DF_Transposition,
    pd,
):
    Oscilation_DF_Baseline['κb;κf'] = 'baseline'
    Oscilation_DF_Transposition['κb;κf'] = '+transposition'
    Oscilation_DF_TPV['κb;κf'] = '+both'

    plotFinal = pd.concat([Oscilation_DF_Baseline, Oscilation_DF_Transposition, Oscilation_DF_TPV])
    return plotFinal,


@app.cell
def __(alt, plasmidN_, plotFinal):
    # First exp comparison, with jumping, vary amp
    t0 = alt.Chart(plotFinal).mark_line(
        point=alt.OverlayMarkDef(size=10, opacity=0.7)).encode(
        x=alt.X('Time', scale=alt.Scale(domain=[-50, 7550]), axis=alt.Axis(labelFontSize=30, values=[0, 2500, 5000, 7500]
                                      , labelAngle=0, ticks=False, grid=False, title=None)),
        y=alt.Y('TCN', scale=alt.Scale(domain=[1, plasmidN_], type='log'), axis=alt.Axis(labelFontSize=30, ticks=False, values=[1,10,plasmidN_], grid=False, title=None)),
        color=alt.Color('κb;κf:O', scale=alt.Scale(domain=['baseline', '+transposition', '+both'], range=['#C3C5DE', '#A2A0C8', '#644F9C']), legend=None)
    ).properties(
        width=400,
        height=150
    ).interactive()

    t0.configure_axis(
        # labelFontSize=14,
        titleFontSize=16
    ).configure_header(
        titleFontSize=20,
        labelFontSize=22 
    ).properties(
        title=alt.TitleParams(
            text= '', #'with selection',
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )
    return t0,


@app.cell
def __(t0):
    t0.save('SF3B_Fluc.pdf')
    return


@app.cell
def __(DownT, UpT, np, pd, plotFinal):
    def calculateLossRateAvg(df):
        return_df_list = []

        for k in list(set(plotFinal['κb;κf'])):
            # print(k)
            currentDF = plotFinal[plotFinal['κb;κf']==k]
            # For baseline case where cells die out after one round of no selection
            if len(currentDF) == DownT+1:

                loss01 = \
                currentDF[(currentDF['Time'] == 0)]['TCN'].values[0] / \
                currentDF[(currentDF['Time'] == DownT)]['TCN'].values[0]

                gen01 = np.log2(currentDF[(currentDF['Time'] == DownT)]['OD'].values[0] / \
                currentDF[(currentDF['Time'] == 0)]['OD'].values[0])

                loss_rate01 = np.log(loss01)/gen01
                loss_rate_List = [loss_rate01]

            # For other cases
            else:
                loss01 = \
                currentDF[(currentDF['Time'] == 0)]['TCN'].values[0] / \
                currentDF[(currentDF['Time'] == DownT)]['TCN'].values[0]               
                gen01 = np.log2(currentDF[(currentDF['Time'] == DownT)]['OD'].values[0] / \
                currentDF[(currentDF['Time'] == 0)]['OD'].values[0])


                loss02 = \
                currentDF[(currentDF['Time'] == DownT+UpT+2)]['TCN'].values[0] / \
                currentDF[(currentDF['Time'] == 2*DownT+UpT+2)]['TCN'].values[0]  
                gen02 = np.log2(currentDF[(currentDF['Time'] == 2*DownT+UpT+2)]['OD'].values[0] / \
                currentDF[(currentDF['Time'] == DownT+UpT+2)]['OD'].values[0])

                loss03 = \
                currentDF[(currentDF['Time'] == 2*DownT+2*UpT+4)]['TCN'].values[0] / \
                currentDF[(currentDF['Time'] == 3*DownT+2*UpT+4)]['TCN'].values[0]          
                gen03 = np.log2(currentDF[(currentDF['Time'] == 3*DownT+2*UpT+4)]['OD'].values[0] / \
                currentDF[(currentDF['Time'] == 2*DownT+2*UpT+4)]['OD'].values[0])

                loss_rate01 = np.log(loss01)/gen01
                loss_rate02 = np.log(loss02)/gen01
                loss_rate03 = np.log(loss03)/gen01

                loss_rate_List = [loss_rate01, loss_rate02, loss_rate03]

            loss_rate_mean = np.mean(loss_rate_List)

            loss_rate_DF = pd.DataFrame({'PCN': [currentDF['PCN'].values[0]],
                                         'κb;κf': [currentDF['κb;κf'].values[0]],
                                         'jumping':[currentDF['jumping'].values[0]],
                                         'plasmidAmp':[currentDF['plasmidAmp'].values[0]],
                                         'lossRate': [loss_rate_mean]})


            return_df_list.append(loss_rate_DF)

            return_df = pd.concat(return_df_list)


        return return_df
    return calculateLossRateAvg,


@app.cell
def __(calculateLossRateAvg, pd, plotFinal):
    LossRatePlot = calculateLossRateAvg(plotFinal)

    # Order the category
    LossRatePlot['κb;κf'] = pd.Series(
                LossRatePlot['κb;κf'], dtype="category")  

    LossRatePlot['κb;κf'] = LossRatePlot[
        'κb;κf'
    ].cat.reorder_categories(
        [
            "baseline",
            "+transposition",
            "+both"
        ],
        ordered=True,       
    )
    return LossRatePlot,


@app.cell
def __(LossRatePlot, alt):
    # loss rate per generation

    Fig2GLeft = alt.Chart(LossRatePlot).mark_point(size=100, filled=True).encode(
        alt.X("κb;κf", axis=alt.Axis(labelFontSize=30, titleFontSize=15, ticks=False, grid=False, labels=False), title= None),
        alt.Y("lossRate", scale=alt.Scale(domain=[-0.01, 0.21]), axis=alt.Axis(labelFontSize=30, tickCount=2, ticks=False, grid=False, title=None)), #scale=alt.Scale(domain=[-0.05, 0.7]), 
        color=alt.Color('κb;κf', scale=alt.Scale(domain=['baseline', '+transposition', '+both'], range=['#C3C5DE', '#A2A0C8', '#644F9C']), legend=None)
    ).properties(
        width=150,
        height=150,
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    ).interactive()

    Fig2GLeft
    return Fig2GLeft,


@app.cell
def __():
    # Fig2GLeft.save('Fig2GLeft.pdf')
    return


@app.cell
def __(DownT, UpT, np, pd, plotFinal):
    def calculateGainRateAvg(df):
        return_df_list = []

        for k in list(set(plotFinal['κb;κf'])):
            # print(k)
            currentDF = plotFinal[plotFinal['κb;κf']==k]

            if len(currentDF) == DownT+1:
                gain_rate_List = [0]

            else:            
                gain01 = \
                currentDF[(currentDF['Time'] == DownT+UpT+1)]['TCN'].values[0] / \
                currentDF[(currentDF['Time'] == DownT+1)]['TCN'].values[0]     
                gen01 = \
                np.log2(currentDF[(currentDF['Time'] == DownT+UpT+1)]['OD'].values[0] / \
                currentDF[(currentDF['Time'] == DownT+1)]['OD'].values[0])

                gain02 = \
                currentDF[(currentDF['Time'] == 2*DownT+2*UpT+3)]['TCN'].values[0] / \
                currentDF[(currentDF['Time'] == 2*DownT+UpT+3)]['TCN'].values[0]          
                gen02 = \
                np.log2(currentDF[(currentDF['Time'] == 2*DownT+2*UpT+3)]['OD'].values[0] / \
                currentDF[(currentDF['Time'] == 2*DownT+1*UpT+3)]['OD'].values[0])

                gain03 = \
                currentDF[(currentDF['Time'] == 3*DownT+3*UpT+5)]['TCN'].values[0] / \
                currentDF[(currentDF['Time'] == 3*DownT+2*UpT+5)]['TCN'].values[0]        
                gen03 = \
                np.log2(currentDF[(currentDF['Time'] == 3*DownT+3*UpT+5)]['OD'].values[0] / \
                currentDF[(currentDF['Time'] == 3*DownT+2*UpT+5)]['OD'].values[0])

                gain_rate01 = np.log(gain01)/gen01
                gain_rate02 = np.log(gain02)/gen02
                gain_rate03 = np.log(gain03)/gen03

                gain_rate_List = [gain_rate01, gain_rate02, gain_rate03]


            gain_rate_mean = np.mean(gain_rate_List)

            gain_rate_DF = pd.DataFrame({'PCN': [currentDF['PCN'].values[0]],
                                         'κb;κf': [currentDF['κb;κf'].values[0]],
                                         'jumping':[currentDF['jumping'].values[0]],
                                         'plasmidAmp':[currentDF['plasmidAmp'].values[0]],
                                         'lossRate': [gain_rate_mean]})

            return_df_list.append(gain_rate_DF)

            return_df = pd.concat(return_df_list)


        return return_df
    return calculateGainRateAvg,


@app.cell
def __(calculateGainRateAvg, pd, plotFinal):
    GainRatePlot = calculateGainRateAvg(plotFinal)
    # Order the category
    GainRatePlot['κb;κf'] = pd.Series(
                GainRatePlot['κb;κf'], dtype="category")  

    GainRatePlot['κb;κf'] = GainRatePlot[
        'κb;κf'
    ].cat.reorder_categories(
        [
            "baseline",
            "+transposition",
            "+both"
        ],
        ordered=True,       
    )
    return GainRatePlot,


@app.cell
def __(GainRatePlot, alt):
    # loss rate per generation

    Fig2GRight = alt.Chart(GainRatePlot).mark_point(size=100, filled=True).encode(
        alt.X("κb;κf", axis=alt.Axis(labelFontSize=30, titleFontSize=15, ticks=False, grid=False, labels=False), title= None),
        alt.Y("lossRate", scale=alt.Scale(domain=[-0.01, 0.11]), axis=alt.Axis(labelFontSize=30, tickCount=1, ticks=False, grid=False, title=None)), #scale=alt.Scale(domain=[-0.05, 0.7]), 
        color=alt.Color('κb;κf', scale=alt.Scale(domain=['baseline', '+transposition', '+both'], range=['#C3C5DE', '#A2A0C8', '#644F9C']), legend=None)
    ).properties(
        width=150,
        height=150,
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).properties(
        title=alt.TitleParams(
            text='',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    ).interactive()

    Fig2GRight
    return Fig2GRight,


@app.cell
def __():
    # Fig2GRight.save('Fig2GRight.pdf')
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
