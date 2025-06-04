import marimo

__generated_with = "0.7.5"
app = marimo.App(width="full")


@app.cell
def __():
    import marimo as mo

    import os
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    from scipy import stats
    import re
    from datetime import datetime
    import itertools
    # plotting packages
    import altair as alt
    return (
        alt,
        datetime,
        itertools,
        matplotlib,
        mo,
        np,
        os,
        pd,
        plt,
        re,
        sns,
        stats,
    )


@app.cell
def __(alt):
    alt.data_transformers.disable_max_rows()
    return


@app.cell
def __():
    figWidth = 600
    figHeight = 270
    colorScheme = 'plasma'

    ## The end of each treatment course
    UP_DAYS = [0, 7, 14, 21]
    DOWN_DAYS = [5, 12, 19]

    TreatmentDilution = 30000
    LossTreatmentDays = 5
    SelectionTreatmentDays = 2
    PlasmidList = ["pBR322", "CloDF", "pUC"] #, "chromosome"
    return (
        DOWN_DAYS,
        LossTreatmentDays,
        PlasmidList,
        SelectionTreatmentDays,
        TreatmentDilution,
        UP_DAYS,
        colorScheme,
        figHeight,
        figWidth,
    )


@app.cell
def __(mo):
    mo.md(r"# Process calibration data first")
    return


@app.cell
def __(pd):
    calibrationFile = pd.ExcelFile('2024-1001_yuanchi.xls')
    return calibrationFile,


@app.cell
def __(calibrationFile):
    calibrationDF = calibrationFile.parse('Results', header = 46).loc[:336, ['Well Position', 'Sample Name', 'Target Name', \
                                                      'Rep', 'Dilution', 'PCN', 'CT']]
    pUC_DF = calibrationDF[calibrationDF['PCN'] == 'pUC']
    p15A_DF = calibrationDF[calibrationDF['PCN'] == 'p15A']

    # Locate entried for each target
    pUC_df_cm = pUC_DF[pUC_DF['Target Name'] == 'cmR']
    pUC_df_kan = pUC_DF[pUC_DF['Target Name'] == 'kanR']
    pUC_df_tet = pUC_DF[pUC_DF['Target Name'] == 'tetA']

    p15A_df_cm = p15A_DF[p15A_DF['Target Name'] == 'cmR']
    p15A_df_kan = p15A_DF[p15A_DF['Target Name'] == 'kanR']
    p15A_df_tet = p15A_DF[p15A_DF['Target Name'] == 'tetA']

    # Use entries with only determined CT values
    pUC_df_cm_used = pUC_df_cm[(pUC_df_cm['Dilution'] != -6) & (pUC_df_cm['Dilution'] != -7)]
    pUC_df_kan_used = pUC_df_kan[(pUC_df_kan['Dilution'] != -6) & (pUC_df_kan['Dilution'] != -7)]
    pUC_df_tet_used = pUC_df_tet[(pUC_df_tet['Dilution'] != -6) & (pUC_df_tet['Dilution'] != -7) ]

    p15A_df_cm_used = p15A_df_cm[p15A_df_cm['Dilution'] != -6]
    p15A_df_kan_used = p15A_df_kan[p15A_df_kan['Dilution'] != -6]
    p15A_df_tet_used = p15A_df_tet[p15A_df_tet['Dilution'] != -6]
    return (
        calibrationDF,
        p15A_DF,
        p15A_df_cm,
        p15A_df_cm_used,
        p15A_df_kan,
        p15A_df_kan_used,
        p15A_df_tet,
        p15A_df_tet_used,
        pUC_DF,
        pUC_df_cm,
        pUC_df_cm_used,
        pUC_df_kan,
        pUC_df_kan_used,
        pUC_df_tet,
        pUC_df_tet_used,
    )


@app.cell
def __(pUC_df_cm_used, pUC_df_kan_used, pUC_df_tet_used, plt, sns):
    figPUC, axesPUC = plt.subplots(1,3, figsize = (7, 2), dpi = 200)
    plt.tight_layout()
    sns.scatterplot(x = 'Dilution', y = 'CT', data = pUC_df_cm_used, ax = axesPUC[0], 
                  marker = 'X')
    sns.scatterplot(x = 'Dilution', y = 'CT', data = pUC_df_tet_used, ax = axesPUC[1], 
                 marker = 'X')
    sns.scatterplot(x = 'Dilution', y = 'CT', data = pUC_df_kan_used, ax = axesPUC[2], 
                 marker = 'X')

    plt.suptitle('pTarget only - standard curves', va = 'bottom')
    axesPUC[0].set_ylabel('cm CT (Chromosome)')
    axesPUC[1].set_ylabel('tet CT (Transposon)')
    axesPUC[2].set_ylabel('kan CT (Plasmid)')
    axesPUC[0].set_xlabel('Log10(Dilution rate)')
    axesPUC[1].set_xlabel('Log10(Dilution rate)')
    axesPUC[2].set_xlabel('Log10(Dilution rate)')
    axesPUC[0].set_xticks([-5,-4,-3,-2,-1,0])
    axesPUC[1].set_xticks([-5,-4,-3,-2,-1,0])
    axesPUC[2].set_xticks([-5,-4,-3,-2,-1,0])
    plt.show()
    return axesPUC, figPUC


@app.cell
def __(p15A_df_cm_used, p15A_df_kan_used, p15A_df_tet_used, plt, sns):
    figP15A, axesP15A = plt.subplots(1,3, figsize = (7, 2), dpi = 200)

    plt.tight_layout()
    sns.scatterplot(x = 'Dilution', y = 'CT', data = p15A_df_cm_used, ax = axesP15A[0], 
                  marker = 'X')
    sns.scatterplot(x = 'Dilution', y = 'CT', data = p15A_df_tet_used, ax = axesP15A[1], 
                 marker = 'X')
    sns.scatterplot(x = 'Dilution', y = 'CT', data = p15A_df_kan_used, ax = axesP15A[2], 
                 marker = 'X')

    plt.suptitle('pTarget only - standard curves (p15A)', va = 'bottom')
    axesP15A[0].set_ylabel('cm CT (Chromosome)')
    axesP15A[1].set_ylabel('tet CT (Transposon)')
    axesP15A[2].set_ylabel('kan CT (Plasmid)')
    axesP15A[0].set_xlabel('Log10(Dilution rate)')
    axesP15A[1].set_xlabel('Log10(Dilution rate)')
    axesP15A[2].set_xlabel('Log10(Dilution rate)')
    axesP15A[0].set_xticks([-5,-4,-3,-2,-1,0])
    axesP15A[1].set_xticks([-5,-4,-3,-2,-1,0])
    axesP15A[2].set_xticks([-5,-4,-3,-2,-1,0])

    plt.show()
    return axesP15A, figP15A


@app.cell
def __(pUC_df_cm_used, pUC_df_kan_used, pUC_df_tet_used, plt, sns):
    plt.figure(figsize=(4, 3))
    sns.scatterplot(x = 'Dilution', y = 'CT', data = pUC_df_cm_used, label='Chromosome', \
                  marker = 'X')
    sns.scatterplot(x = 'Dilution', y = 'CT', data = pUC_df_tet_used, label='Transposon', \
                 marker = 'X')
    sns.scatterplot(x = 'Dilution', y = 'CT', data = pUC_df_kan_used, label='Plasmid', \
                 marker = 'X')

    plt.xlabel('Log10(Dilution rate)')
    plt.ylabel('CT')
    plt.title('pTarget only - standard curves (pUC)')
    plt.show()
    return


@app.cell
def __(p15A_df_cm_used, p15A_df_kan_used, p15A_df_tet_used, plt, sns):
    plt.figure(figsize=(4, 3))
    sns.scatterplot(x = 'Dilution', y = 'CT', data = p15A_df_cm_used, label='Chromosome', \
                  marker = 'X')
    sns.scatterplot(x = 'Dilution', y = 'CT', data = p15A_df_tet_used, label='Transposon', \
                 marker = 'X')
    sns.scatterplot(x = 'Dilution', y = 'CT', data = p15A_df_kan_used, label='Plasmid', \
                 marker = 'X')

    plt.xlabel('Log10(Dilution rate)')
    plt.ylabel('CT')
    plt.title('pTarget only - standard curves (p15A)')
    plt.show()
    return


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"# Long-term experimental data analysis starting here")
    return


@app.cell
def __(np):
    def calculateSlopeAmpEff(targetedDFList):
        slopeList = []
        amplificationFactorList = []
        amplificationEfficiencyList = []

        for df in targetedDFList:
            slope = np.polynomial.polynomial.polyfit(list(df['Dilution']), list(df['CT']), 1)[1]
            amplificationFactor = 10**(-1*(1/slope))
            amplificationEfficiency = (amplificationFactor-1)*100

            slopeList.append(slope)
            amplificationFactorList.append(amplificationFactor)
            amplificationEfficiencyList.append(amplificationEfficiency)

        return slopeList, amplificationFactorList, amplificationEfficiencyList
    return calculateSlopeAmpEff,


@app.cell
def __(
    calculateSlopeAmpEff,
    pUC_df_cm_used,
    pUC_df_kan_used,
    pUC_df_tet_used,
):
    slopeList, ampFactorList, effList = calculateSlopeAmpEff([pUC_df_cm_used, pUC_df_kan_used, pUC_df_tet_used])
    return ampFactorList, effList, slopeList


@app.cell
def __(mo):
    mo.md(r"# Import daily qPCR results")
    return


@app.cell
def __(pd):
    D0 = pd.ExcelFile('D0-2024-10-05.xls').parse('Results', header = 46).loc[:588, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D0_noKan = pd.ExcelFile('D0_noKan-2024-11-08.xls').parse('Results', header = 46).loc[:288, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D1 = pd.ExcelFile('D1-2024-10-15.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D1_noKan = pd.ExcelFile('D1_noKan-2024-10-16.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D2 = pd.ExcelFile('D2-2024-10-25_02.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D2_noKan = pd.ExcelFile('D2_noKan-2024-10-26.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D3 = pd.ExcelFile('D3-2024-10-28.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D3_noKan = pd.ExcelFile('D3_noKan-2024-10-28.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D4 = pd.ExcelFile('D4-2024-11-02.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D4_noKan = pd.ExcelFile('D4_noKan-2024-11-03.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D5 = pd.ExcelFile('D5-2024-11-25.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D5_noKan = pd.ExcelFile('D5_noKan-2024-10-19.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D6 = pd.ExcelFile('D6-2024-10-30.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D6_noKan = pd.ExcelFile('D6_noKan-2024-10-30.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D7 = pd.ExcelFile('D7-2024-10-29_02.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D7_noKan = pd.ExcelFile('D7_noKan-2024-10-29.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D8 = pd.ExcelFile('D8-2024-11-03.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D8_noKan = pd.ExcelFile('D8_noKan-2024-11-03.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D9 = pd.ExcelFile('D9-2024-11-06.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D9_noKan = pd.ExcelFile('D9_noKan-2024-11-06.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D10 = pd.ExcelFile('D10-2024-11-07.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D10_noKan = pd.ExcelFile('D10_noKan-2024-11-07.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D11 = pd.ExcelFile('D11-2024-11-07.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D11_noKan = pd.ExcelFile('D11_noKan-2024-11-07.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D12 = pd.ExcelFile('D12-2024-11-10.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D12_noKan = pd.ExcelFile('D12_noKan-2024-11-08.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D13 = pd.ExcelFile('D13-2024-11-08.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D13_noKan = pd.ExcelFile('D13_noKan-2024-11-08.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D14 = pd.ExcelFile('D14-2024-10-13.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D14_noKan = pd.ExcelFile('D14_noKan-2024-10-17.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D15 = pd.ExcelFile('D15-2024-11-08.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D15_noKan = pd.ExcelFile('D15_noKan-2024-11-08.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D16 = pd.ExcelFile('D16-2024-11-06.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D16_noKan = pd.ExcelFile('D16_noKan-2024-11-06.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D17 = pd.ExcelFile('D17-2024-11-10.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D17_noKan = pd.ExcelFile('D17_noKan-2024-11-10.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D18 = pd.ExcelFile('D18-2024-11-04.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D18_noKan = pd.ExcelFile('D18_noKan-2024-11-04_02.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D19 = pd.ExcelFile('D19-2024-10-11.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D19_noKan = pd.ExcelFile('D19_noKan-2024-10-17.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D20 = pd.ExcelFile('D20-2024-10-14.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D20_noKan = pd.ExcelFile('D20_noKan-2024-10-16.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]

    D21 = pd.ExcelFile('D21-2024-10-30.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    D21_noKan = pd.ExcelFile('D21_noKan-2024-10-30.xls').parse('Results', header = 46).loc[:336, ['Well Position', \
                                                      'qPCRRep', 'ExpRep', 'PCN', 'Target Name', 'CT']]
    return (
        D0,
        D0_noKan,
        D1,
        D10,
        D10_noKan,
        D11,
        D11_noKan,
        D12,
        D12_noKan,
        D13,
        D13_noKan,
        D14,
        D14_noKan,
        D15,
        D15_noKan,
        D16,
        D16_noKan,
        D17,
        D17_noKan,
        D18,
        D18_noKan,
        D19,
        D19_noKan,
        D1_noKan,
        D2,
        D20,
        D20_noKan,
        D21,
        D21_noKan,
        D2_noKan,
        D3,
        D3_noKan,
        D4,
        D4_noKan,
        D5,
        D5_noKan,
        D6,
        D6_noKan,
        D7,
        D7_noKan,
        D8,
        D8_noKan,
        D9,
        D9_noKan,
    )


@app.cell
def __(
    D0,
    D0_noKan,
    D1,
    D10,
    D10_noKan,
    D11,
    D11_noKan,
    D12,
    D12_noKan,
    D13,
    D13_noKan,
    D14,
    D14_noKan,
    D15,
    D15_noKan,
    D16,
    D16_noKan,
    D17,
    D17_noKan,
    D18,
    D18_noKan,
    D19,
    D19_noKan,
    D1_noKan,
    D2,
    D20,
    D20_noKan,
    D21,
    D21_noKan,
    D2_noKan,
    D3,
    D3_noKan,
    D4,
    D4_noKan,
    D5,
    D5_noKan,
    D6,
    D6_noKan,
    D7,
    D7_noKan,
    D8,
    D8_noKan,
    D9,
    D9_noKan,
):
    datafiles_Kan = [(D0,0), (D1,1), (D2,2), (D3,3), (D4,4), (D5,5), (D6,6), (D7,7), (D8,8), (D9,9), (D10,10),\
                     (D11,11), (D12,12), (D13,13), (D14,14), (D15,15), (D16,16), (D17,17), (D18,18), \
                     (D19,19), (D20,20), (D21,21)]

    datafiles_noKan = [(D0_noKan,0), (D1_noKan,1), (D2_noKan,2), (D3_noKan,3), (D4_noKan,4), (D5_noKan,5), (D6_noKan,6), \
                       (D7_noKan,7), (D8_noKan,8), (D9_noKan,9), (D10_noKan,10), (D11_noKan,11), (D12_noKan,12), \
                       (D13_noKan,13), (D14_noKan,14), (D15_noKan,15), (D16_noKan,16), (D17_noKan,17), \
                       (D18_noKan,18), (D19_noKan,19), (D20_noKan,20), (D21_noKan,21)]
    return datafiles_Kan, datafiles_noKan


@app.cell
def __(ampFactorList, effList, itertools, pd, slopeList):
    def processOneDay(DF, day):
        print(day)
        calculatedDF = pd.DataFrame()

        # Get all calibrated parameters
        DF.loc[DF['Target Name'] == 'cmR', 'Slope'] = slopeList[0]
        DF.loc[DF['Target Name'] == 'kanR', 'Slope'] = slopeList[1]
        DF.loc[DF['Target Name'] == 'tetA', 'Slope'] = slopeList[2]

        DF.loc[DF['Target Name'] == 'cmR', 'AmpFactor'] = ampFactorList[0]
        DF.loc[DF['Target Name'] == 'kanR', 'AmpFactor'] = ampFactorList[1]
        DF.loc[DF['Target Name'] == 'tetA', 'AmpFactor'] = ampFactorList[2]

        DF.loc[DF['Target Name'] == 'cmR', 'Eff'] = effList[0]
        DF.loc[DF['Target Name'] == 'kanR', 'Eff'] = effList[1]
        DF.loc[DF['Target Name'] == 'tetA', 'Eff'] = effList[2]

        qPCRRepList = list(set(DF['qPCRRep']))
        ExpRepList = list(set(DF['ExpRep']))
        PCNList = list(set(DF['PCN']))
        PCNList = [i for i in PCNList if (i != 'H2O') & (i != 'Blank')]

        qPCR_Exp_PCN_Tuple_List = list(itertools.product(qPCRRepList, ExpRepList, PCNList))

        for t in qPCR_Exp_PCN_Tuple_List:
            df_current = DF[(DF['qPCRRep'] == t[0]) & (DF['ExpRep'] == t[1]) & (DF['PCN'] == t[2])]
            currentSlice = DF.loc[df_current.index]

            if currentSlice.empty == False:
                if 'tetA' in (list(currentSlice['Target Name'])):
                    currentSlice = currentSlice[currentSlice['CT'] != 'Undetermined']

                    cmR_CT = pd.DataFrame(currentSlice[currentSlice['Target Name'] == 'cmR']['CT']).values[0][0]
                    tetA_CT = pd.DataFrame(currentSlice[currentSlice['Target Name'] == 'tetA']['CT']).values[0][0]

                    cmR_ampFactor = currentSlice[currentSlice['Target Name'] == 'cmR']['AmpFactor'].values[0]
                    tetA_ampFactor = currentSlice[currentSlice['Target Name'] == 'tetA']['AmpFactor'].values[0]

                    TCopy = cmR_ampFactor**(cmR_CT)/tetA_ampFactor**(tetA_CT)
                    currentSlice['transposon copy'] = TCopy

                    if currentSlice['PCN'].values[0] != 'Chromosome':
                        kanR_CT = pd.DataFrame(currentSlice[currentSlice['Target Name'] == 'kanR']['CT']).values[0][0]
                        kanR_ampFactor = currentSlice[currentSlice['Target Name'] == 'kanR']['AmpFactor'].values[0]
                        PCopy = cmR_ampFactor**(cmR_CT)/kanR_ampFactor**(kanR_CT)
                        currentSlice['plasmid copy'] = PCopy
                        currentSlice['T/P'] = (TCopy-1)/PCopy
                    else:
                        currentSlice['plasmid copy'] = 0

                else:
                    currentSlice = currentSlice[currentSlice['CT'] != 'Undetermined']

                    cmR_CT = pd.DataFrame(currentSlice[currentSlice['Target Name'] == 'cmR']['CT']).values[0][0]
                    kanR_CT = pd.DataFrame(currentSlice[currentSlice['Target Name'] == 'kanR']['CT']).values[0][0]

                    cmR_ampFactor = currentSlice[currentSlice['Target Name'] == 'cmR']['AmpFactor'].values[0]
                    kanR_ampFactor = currentSlice[currentSlice['Target Name'] == 'kanR']['AmpFactor'].values[0]

                    currentSlice['transposon copy'] = 0
                    currentSlice['T/P'] = 0

                    if kanR_CT != 'Undetermined':
                        currentSlice['plasmid copy'] = cmR_ampFactor**(cmR_CT)/kanR_ampFactor**(kanR_CT)
                    else:
                        currentSlice['plasmid copy'] = 0

                calculatedDF = pd.concat([calculatedDF, currentSlice], ignore_index=True)

        calculatedDF["Day"] = day
        calculatedDF['qPCRRep'] = calculatedDF['qPCRRep'].astype(int)  
        calculatedDF['ExpRep'] = calculatedDF['ExpRep'].astype(int).astype(str)


        return calculatedDF
    return processOneDay,


@app.cell
def __(datafiles_Kan, datafiles_noKan, pd, processOneDay):
    dataframes_tet10_kan = \
    [processOneDay(pair[0], pair[1]) for pair in datafiles_Kan]
    dataframes_tet10_noKan = \
    [processOneDay(pair[0], pair[1]) for pair in datafiles_noKan]

    completeDF_tet10_kan = pd.concat(dataframes_tet10_kan, ignore_index=True)
    completeDF_tet10_kan['tet'] = 10
    completeDF_tet10_noKan = pd.concat(dataframes_tet10_noKan, ignore_index=True)
    completeDF_tet10_noKan['tet'] = 10
    return (
        completeDF_tet10_kan,
        completeDF_tet10_noKan,
        dataframes_tet10_kan,
        dataframes_tet10_noKan,
    )


@app.cell
def __(completeDF_tet10_kan, completeDF_tet10_noKan, pd):
    # Take out all the calibration and testing wells for further analysis
    tet10_kan_DF = completeDF_tet10_kan[(completeDF_tet10_kan['PCN'] != 'cyro_chromosome') & \
                                   (completeDF_tet10_kan['PCN'] != 'cyro_pSC101') & \
                                   (completeDF_tet10_kan['PCN'] != 'p15ASep')]

    tet10_kan_DF["PCN"] = pd.Series(
        tet10_kan_DF["PCN"], dtype="category"
    )

    tet10_kan_DF["PCN"] = tet10_kan_DF[
                "PCN"
            ].cat.reorder_categories(
                [
                    "p15A",
                    "pBR322",
                    "CloDF",
                    "pUC",
                    "Chromosome",
                ],
                ordered=True,
            )


    tet10_noKan_DF = completeDF_tet10_noKan[(completeDF_tet10_noKan['PCN'] != 'cyro_chromosome') & \
                                   (completeDF_tet10_noKan['PCN'] != 'cyro_pSC101') & \
                                   (completeDF_tet10_noKan['PCN'] != 'p15ASep')]

    tet10_noKan_DF["PCN"] = pd.Series(
        tet10_noKan_DF["PCN"], dtype="category"
    )

    tet10_noKan_DF["PCN"] = tet10_noKan_DF[
                "PCN"
            ].cat.reorder_categories(
                [
                    "p15A",
                    "pBR322",
                    "CloDF",
                    "pUC",
                ],
                ordered=True,
            )
    return tet10_kan_DF, tet10_noKan_DF


@app.cell
def __(np, tet10_kan_DF, tet10_noKan_DF):
    tet10_kan_DF['log(plasmid copy)'] = np.log10(tet10_kan_DF['plasmid copy'])
    tet10_kan_DF['log(transposon copy)'] = np.log10(tet10_kan_DF['transposon copy'])

    tet10_noKan_DF['log(plasmid copy)'] = np.log10(tet10_noKan_DF['plasmid copy'])
    tet10_noKan_DF['log(transposon copy)'] = np.log10(tet10_noKan_DF['transposon copy'])
    return


@app.cell
def __(tet10_kan_DF, tet10_noKan_DF):
    completeDF_tet10_kan_copy = tet10_kan_DF.copy(deep=True)
    completeDF_tet10_noKan_copy = tet10_noKan_DF.copy(deep=True)
    return completeDF_tet10_kan_copy, completeDF_tet10_noKan_copy


@app.cell
def __(tet10_kan_DF, tet10_noKan_DF):
    tet10_kan_DF_Use = \
    tet10_kan_DF[(tet10_kan_DF['PCN'] != 'p15A')] # & (tet10_kan_DF['PCN'] != 'Chromosome')

    tet10_noKan_DF_Use = \
    tet10_noKan_DF[(tet10_noKan_DF['PCN'] != 'p15A')]

    tet10_kan_DF_Use_pUC = tet10_kan_DF_Use[tet10_kan_DF_Use['PCN'] == 'pUC']
    tet10_noKan_DF_Use_pUC = tet10_noKan_DF_Use[tet10_noKan_DF_Use['PCN'] == 'pUC']
    return (
        tet10_kan_DF_Use,
        tet10_kan_DF_Use_pUC,
        tet10_noKan_DF_Use,
        tet10_noKan_DF_Use_pUC,
    )


@app.cell
def __(tet10_kan_DF_Use, tet10_noKan_DF_Use):
    tet10_kan_DF_Use_grouped = tet10_kan_DF_Use.groupby(['Day', 'ExpRep', 'PCN'])[['log(transposon copy)', 'log(plasmid copy)', 'transposon copy', 'plasmid copy', 'T/P']].mean().reset_index().dropna()

    tet10_noKan_DF_Use_grouped = tet10_noKan_DF_Use.groupby(['Day', 'ExpRep', 'PCN'])[['log(transposon copy)', 'log(plasmid copy)', 'transposon copy', 'plasmid copy', 'T/P']].mean().reset_index().dropna()
    return tet10_kan_DF_Use_grouped, tet10_noKan_DF_Use_grouped


@app.cell
def __(pd):
    # Get the control chromosome strain results from another file
    noTransposon_kan_qPCR_timeCourse_DF = pd.read_csv('noTransposon_kan_qPCR_timeCourse.csv')
    noTransposon_kan_qPCR_timeCourse_DF_Chr = noTransposon_kan_qPCR_timeCourse_DF[noTransposon_kan_qPCR_timeCourse_DF['PCN'] == 'Chromosome']

    noTransposon_noKan_qPCR_timeCourse_DF = pd.read_csv('noTransposon_noKan_qPCR_timeCourse.csv')
    noTransposon_noKan_qPCR_timeCourse_DF_Chr = noTransposon_noKan_qPCR_timeCourse_DF[noTransposon_noKan_qPCR_timeCourse_DF['PCN'] == 'Chromosome']
    return (
        noTransposon_kan_qPCR_timeCourse_DF,
        noTransposon_kan_qPCR_timeCourse_DF_Chr,
        noTransposon_noKan_qPCR_timeCourse_DF,
        noTransposon_noKan_qPCR_timeCourse_DF_Chr,
    )


@app.cell
def __(
    noTransposon_kan_qPCR_timeCourse_DF_Chr,
    noTransposon_noKan_qPCR_timeCourse_DF_Chr,
    pd,
    tet10_kan_DF_Use_grouped,
    tet10_noKan_DF_Use_grouped,
):
    # For the control
    tet10_kan_DF_Use_grouped_final = pd.concat([tet10_kan_DF_Use_grouped, noTransposon_kan_qPCR_timeCourse_DF_Chr], ignore_index=True)
    tet10_noKan_DF_Use_grouped_final = pd.concat([tet10_noKan_DF_Use_grouped, noTransposon_noKan_qPCR_timeCourse_DF_Chr], ignore_index=True)
    return tet10_kan_DF_Use_grouped_final, tet10_noKan_DF_Use_grouped_final


@app.cell
def __(
    pd,
    tet10_kan_DF_Use_grouped_final,
    tet10_noKan_DF_Use_grouped_final,
):
    # Order dataframes
    tet10_kan_DF_Use_grouped_final["PCN"] = pd.Series(
            tet10_kan_DF_Use_grouped_final["PCN"], dtype="category"
        )

    tet10_kan_DF_Use_grouped_final["PCN"] = tet10_kan_DF_Use_grouped_final[
                "PCN"
            ].cat.reorder_categories(
                [
                    "pBR322",
                    "CloDF",
                    "pUC",
                    "Chromosome",
                ],
                ordered=True,
            )

    tet10_noKan_DF_Use_grouped_final["PCN"] = pd.Series(
            tet10_noKan_DF_Use_grouped_final["PCN"], dtype="category"
        )

    tet10_noKan_DF_Use_grouped_final["PCN"] = tet10_noKan_DF_Use_grouped_final[
                "PCN"
            ].cat.reorder_categories(
                [
                    "pBR322",
                    "CloDF",
                    "pUC",
                    "Chromosome",
                ],
                ordered=True,
            )
    return


@app.cell
def __(plt, sns, tet10_kan_DF_Use_grouped):
    ## Experimental conditions: +Kan

    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    palette02 = ["#7316A2", "#EA9953", "#BC5078", ]

    # Create line plot with points
    ax_kan = sns.lineplot(
        data=tet10_kan_DF_Use_grouped,
        x='Day',
        y='log(plasmid copy)',
        hue='PCN',  
        units='ExpRep',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=0.7,       # Match opacity from Altair
        palette=palette02,
        legend=False     # No legend as in original
    )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(-0.2, 4.5)

    # Configure axis appearance
    ax_kan.tick_params(axis='x', labelsize=25, length=0)
    ax_kan.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([0, 2, 4])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('ln(PCN)', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax_kan, palette02


@app.cell
def __(palette02, plt, sns, tet10_kan_DF_Use_grouped):
    ## Experimental conditions: +Kan

    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    # Create line plot with points
    ax_kan_TCN = sns.lineplot(
        data=tet10_kan_DF_Use_grouped,
        x='Day',
        y='log(transposon copy)',
        hue='PCN',  
        units='ExpRep',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=0.7,       # Match opacity from Altair
        palette=palette02,
        legend=False     # No legend as in original
    )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(-0.2, 4.5)

    # Configure axis appearance
    ax_kan_TCN.tick_params(axis='x', labelsize=25, length=0)
    ax_kan_TCN.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([0, 2, 4])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('ln(TCN)', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax_kan_TCN,


@app.cell
def __(palette02, plt, sns, tet10_kan_DF_Use_grouped):
    ## Experimental conditions: +Kan

    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    # Create line plot with points
    ax_kan_TP = sns.lineplot(
        data=tet10_kan_DF_Use_grouped,
        x='Day',
        y='T/P',
        hue='PCN',  
        units='ExpRep',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=0.7,       # Match opacity from Altair
        palette=palette02,
        legend=False     # No legend as in original
    )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(-0.1, 1.5)

    # Configure axis appearance
    ax_kan_TP.tick_params(axis='x', labelsize=25, length=0)
    ax_kan_TP.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([0, 1])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('T/P', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax_kan_TP,


@app.cell
def __():
    return


@app.cell
def __(palette02, plt, sns, tet10_noKan_DF_Use_grouped):
    ## Experimental conditions: -Kan

    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    # Create line plot with points
    ax_noKan = sns.lineplot(
        data=tet10_noKan_DF_Use_grouped,
        x='Day',
        y='log(plasmid copy)',
        hue='PCN',  
        units='ExpRep',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=0.7,       # Match opacity from Altair
        palette=palette02,
        legend=False     # No legend as in original
    )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(-0.2, 4.5)

    # Configure axis appearance
    ax_noKan.tick_params(axis='x', labelsize=25, length=0)
    ax_noKan.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([0, 2, 4])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('ln(PCN)', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax_noKan,


@app.cell
def __(palette02, plt, sns, tet10_noKan_DF_Use_grouped):
    ## Experimental conditions: -Kan

    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    # Create line plot with points
    ax_noKan_TCN = sns.lineplot(
        data=tet10_noKan_DF_Use_grouped,
        x='Day',
        y='log(transposon copy)',
        hue='PCN',  
        units='ExpRep',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=0.7,       # Match opacity from Altair
        palette=palette02,
        legend=False     # No legend as in original
    )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(-0.2, 4.5)

    # Configure axis appearance
    ax_noKan_TCN.tick_params(axis='x', labelsize=25, length=0)
    ax_noKan_TCN.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([0, 2, 4])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('ln(TCN)', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax_noKan_TCN,


@app.cell
def __(palette02, plt, sns, tet10_noKan_DF_Use_grouped):
    ## Experimental conditions: -Kan

    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    # Create line plot with points
    ax_noKan_TP = sns.lineplot(
        data=tet10_noKan_DF_Use_grouped,
        x='Day',
        y='T/P',
        hue='PCN',  
        units='ExpRep',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=0.7,       # Match opacity from Altair
        palette=palette02,
        legend=False     # No legend as in original
    )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(-0.1, 1.5)

    # Configure axis appearance
    ax_noKan_TP.tick_params(axis='x', labelsize=25, length=0)
    ax_noKan_TP.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([0, 1])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('T/P', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax_noKan_TP,


@app.cell
def __(plt, sns, tet10_kan_DF_Use_grouped_final):
    ## Experimental conditions: +Kan & + control group

    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    palette = ["#EA9953", "#BC5078", "#7316A2", "#D5D5D5"]

    # Create line plot with points
    ax = sns.lineplot(
        data=tet10_kan_DF_Use_grouped_final,
        x='Day',
        y='log(transposon copy)',
        hue='PCN',  
        units='ExpRep',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=1,       # Match opacity from Altair
        palette=palette ,
        legend=False     # No legend as in original
    )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(-5, 5)

    # Configure axis appearance
    ax.tick_params(axis='x', labelsize=25, length=0)
    ax.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([-5, 0, 5])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('ln(TCN)', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax, palette


@app.cell
def __(pd, tet10_kan_DF_Use_grouped, tet10_noKan_DF_Use_grouped):
    # Calculate PCN by averaging over everyday
    tet10_kan_DF_Use_grouped.groupby(['PCN'])['plasmid copy'].mean()

    # Combine 2 dataframes for further analysis
    tet10_kan_DF_Use_grouped['Kan'] = 'Kan'
    tet10_noKan_DF_Use_grouped['Kan'] = 'noKan'

    combined_data = pd.concat([tet10_kan_DF_Use_grouped, tet10_noKan_DF_Use_grouped])
    return combined_data,


@app.cell
def __(mo):
    mo.md(r"# Calculate loss rates")
    return


@app.cell
def __(
    DOWN_DAYS,
    LossTreatmentDays,
    PlasmidList,
    TreatmentDilution,
    UP_DAYS,
    np,
):
    def calculateLossRate(df, qPCRtarget):
        # Make pairs of (Up day, Down day)
        up_down_pair = list(zip(UP_DAYS[:], DOWN_DAYS[:]))
        # print(up_down_pair)

        # Process each plasmid
        for plasmid in PlasmidList:
            if plasmid != 'p15A':
                # Process each pair
                for pair in up_down_pair:
                    # Selection Day result
                    denominator = df[(df['Day'] == pair[0]) & (df['PCN'] == plasmid)][qPCRtarget]
                    # The last day of no selection After the selection day
                    numerator = df[(df['Day'] == pair[1]) & (df['PCN'] == plasmid)][qPCRtarget]
                    
                    # Ensure the indexes align
                    denominator = denominator.reset_index(drop=True)
                    numerator = numerator.reset_index(drop=True)

                    # To make the value positive
                    loss = denominator/numerator

                    # print(loss)

                    loss = loss.dropna()

                    # print(loss)
                    indices = df[(df['Day'] == pair[1]) & (df['PCN'] == plasmid)].index
                    # print(len(indices), len(loss.values))

                    # Assign the calculated loss back
                    df.loc[indices, 'loss'] = loss.values

                    # Take mean of all days  
                    df.loc[indices, 'mean of loss'] = np.mean(df.loc[indices, 'loss'])

                    # Calculate the rate & stats
                    generation = np.log2(TreatmentDilution^LossTreatmentDays)    
                    loss_rate = np.log(loss)/generation
                    loss_rate_mean = np.mean(loss_rate)
                    loss_rate_std = np.std(loss_rate)

                    # Assign calculated loss rate & stats values to the df
                    df.loc[indices, 'rate_of_loss_log('+qPCRtarget+')'] = loss_rate.values
                    df.loc[indices, 'mean_of_loss_rate(log('+qPCRtarget+'))'] = loss_rate_mean
                    df.loc[indices, 'std_of_loss_rate(log('+qPCRtarget+'))'] = loss_rate_std
                    df.loc[indices, 'ci_lower('+qPCRtarget+')'] = loss_rate_mean - loss_rate_std
                    df.loc[indices, 'ci_upper('+qPCRtarget+')'] = loss_rate_mean + loss_rate_std

                    # print(df)

        # df['PCN'] = pd.Series(df['PCN'], dtype="category")
        # # print(set(df['PCN']))

        # df["PCN"] = df[
        #         "PCN"
        #     ].cat.reorder_categories(
        #         [
        #             'pBR322',
        #             'CloDF',
        #             'pUC',
        #         ],
        #         ordered=True,
        #     )

        df = df.fillna(np.nan)

        return df
    return calculateLossRate,


@app.cell
def __(calculateLossRate, pd):
    def makeLossPlotDFs(currentCompleteDF):
        currentCompleteDF_lossCopy = currentCompleteDF.copy(deep=True)
        # for the transposon result
        currentCompleteDF_loss_transposon_df = calculateLossRate(currentCompleteDF_lossCopy, 'transposon copy')
        # for the plasmid result
        currentCompleteDF_loss_plasmid_df = calculateLossRate(currentCompleteDF_lossCopy, 'plasmid copy')
        # for the T/P result
        currentCompleteDF_loss_ratio_df = calculateLossRate(currentCompleteDF_lossCopy, 'T/P')

        # To postprocess for ploting
        currentCompleteDF_loss_transposon_df_cleaned_scatter = currentCompleteDF_loss_transposon_df.groupby(['PCN','ExpRep'])[['rate_of_loss_log(transposon copy)']].mean().reset_index().dropna()

        currentCompleteDF_loss_plasmid_df_cleaned_scatter = currentCompleteDF_loss_plasmid_df.groupby(['PCN','ExpRep'])[['rate_of_loss_log(plasmid copy)']].mean().reset_index().dropna()

        currentCompleteDF_loss_ratio_df_cleaned_scatter = currentCompleteDF_loss_ratio_df.groupby(['PCN','ExpRep'])[['rate_of_loss_log(T/P)']].mean().reset_index().dropna()

        # Identify overlapping column names
        overlapping_columns = set(currentCompleteDF_loss_transposon_df_cleaned_scatter.columns) & set(currentCompleteDF_loss_plasmid_df_cleaned_scatter.columns) & set(currentCompleteDF_loss_ratio_df_cleaned_scatter.columns)

        # Concatenate DataFrames horizontally, keeping only one occurrence of overlapping columns
        result_df = pd.concat([currentCompleteDF_loss_transposon_df_cleaned_scatter, currentCompleteDF_loss_plasmid_df_cleaned_scatter.loc[:, ~currentCompleteDF_loss_plasmid_df_cleaned_scatter.columns.isin(overlapping_columns)]], axis=1)

        result_df = pd.concat([result_df, currentCompleteDF_loss_ratio_df_cleaned_scatter.loc[:, ~currentCompleteDF_loss_ratio_df_cleaned_scatter.columns.isin(overlapping_columns)]], axis=1)

        return result_df
    return makeLossPlotDFs,


@app.cell
def __(
    makeLossPlotDFs,
    pd,
    tet10_kan_DF_Use_grouped,
    tet10_noKan_DF_Use_grouped,
):
    scatter_kan = makeLossPlotDFs(tet10_kan_DF_Use_grouped)
    scatter_noKan = makeLossPlotDFs(tet10_noKan_DF_Use_grouped)

    scatter_kan['Kan'] = 1
    scatter_noKan['Kan'] = 0

    combined_scatter_DF = pd.concat([scatter_kan, scatter_noKan])

    combined_scatter_DF['PCN'] = pd.Series(combined_scatter_DF['PCN'], dtype="category")
    combined_scatter_DF["PCN"] = combined_scatter_DF["PCN"].cat.reorder_categories(
                    [
                    'pBR322',
                    'CloDF',
                    'pUC',
                    ],
                ordered=True,
                    )
    return combined_scatter_DF, scatter_kan, scatter_noKan


@app.cell
def __(alt, colorScheme, combined_scatter_DF):
    # Transposon loss rate per generation

    alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title=None),
        alt.Y("rate_of_loss_log(transposon copy)", scale=alt.Scale(domain=[-0.02, 0.33]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, ticks=False, grid=False), title=None),
        color=alt.Color('PCN', ).scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        width=150,
        height=150
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).properties(
        title=alt.TitleParams(
            text='Transposon loss rate',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    ).interactive()

    return


@app.cell
def __(alt, colorScheme, combined_scatter_DF):
    # plasmid loss rate per generation, PCN

    alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title=None),
        alt.Y("rate_of_loss_log(plasmid copy)", axis=alt.Axis(labelFontSize=30, titleFontSize=15, ticks=False, grid=False), title=None),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        width=150, 
        height=150
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).properties(
        title=alt.TitleParams(
            text='Plasmid loss rate',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    ).interactive()

    # .facet(
    #     column=alt.Column('Kan', header=alt.Header(title=None))
    # )
    return


@app.cell
def __(alt, colorScheme, combined_scatter_DF):
    # T/P loss rate per generation

    alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title=None),
        alt.Y("rate_of_loss_log(T/P)", scale=alt.Scale(type='log'), axis=alt.Axis(labelFontSize=30, titleFontSize=15, ticks=False, grid=False), title=None),

        color=alt.Color('PCN',).scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        width=150,
        height=150
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).properties(
        title=alt.TitleParams(
            text='T/P loss rate',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    ).interactive()

    # .facet(
    #     column=alt.Column('Kan', header=alt.Header(title=None))
    # )
    return


@app.cell
def __(combined_scatter_DF, pd, tet10_kan_DF_Use_grouped):
    averagePCN = pd.DataFrame(tet10_kan_DF_Use_grouped.groupby(['PCN'])['plasmid copy'].mean()).reset_index()

    combined_scatter_DF.loc[combined_scatter_DF['PCN'] == 'pBR322', 'averagePCN'] = averagePCN[averagePCN['PCN'] == 'pBR322']['plasmid copy'].values[0]
    combined_scatter_DF.loc[combined_scatter_DF['PCN'] == 'CloDF', 'averagePCN'] = averagePCN[averagePCN['PCN'] == 'CloDF']['plasmid copy'].values[0]
    combined_scatter_DF.loc[combined_scatter_DF['PCN'] == 'pUC', 'averagePCN'] = averagePCN[averagePCN['PCN'] == 'pUC']['plasmid copy'].values[0]
    return averagePCN,


@app.cell
def __(alt, colorScheme, combined_scatter_DF):
    TCN_loss = alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_loss_log(transposon copy)", scale=alt.Scale(domain=[0, 0.35]), axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='loss rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='TCN',
        width=200,
        height=200
    )

    PCN_loss = alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=25, titleFontSize=25, grid=False), title='PCN'),
        alt.Y("rate_of_loss_log(plasmid copy)", scale=alt.Scale(domain=[0, 0.35]), axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='loss rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='PCN',
        width=200,
        height=200
    )

    Ratio_loss = alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_loss_log(T/P)", scale=alt.Scale(domain=[0, 0.35]), axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='loss rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='T/P',
        width=200,
        height=200
    )

    Loss_combined_chart = TCN_loss | PCN_loss | Ratio_loss

    Loss_combined_chart.interactive()
    return Loss_combined_chart, PCN_loss, Ratio_loss, TCN_loss


@app.cell
def __(alt, colorScheme, combined_scatter_DF):
    TCN_loss_PCNasX = alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("averagePCN", axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_loss_log(transposon copy)", scale=alt.Scale(domain=[0, 0.33]), axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='loss rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='TCN',
        width=200,
        height=200
    )

    PCN_loss_PCNasX = alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("averagePCN", axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_loss_log(plasmid copy)", scale=alt.Scale(domain=[0, 0.33]), axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='loss rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='PCN',
        width=200,
        height=200
    )

    Ratio_loss_PCNasX = alt.Chart(combined_scatter_DF).mark_point(size=70).encode(
        alt.X("averagePCN", axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_loss_log(T/P)", scale=alt.Scale(domain=[0, 0.33]), axis=alt.Axis(labelFontSize=25, titleFontSize=15, grid=False), title='loss rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='T/P',
        width=200,
        height=200
    )

    Loss_combined_PCNasX_chart = TCN_loss_PCNasX | PCN_loss_PCNasX | Ratio_loss_PCNasX

    Loss_combined_PCNasX_chart.interactive()
    return (
        Loss_combined_PCNasX_chart,
        PCN_loss_PCNasX,
        Ratio_loss_PCNasX,
        TCN_loss_PCNasX,
    )


@app.cell
def __(mo):
    mo.md(r"# Calculate gain rates")
    return


@app.cell
def __(
    DOWN_DAYS,
    PlasmidList,
    SelectionTreatmentDays,
    TreatmentDilution,
    UP_DAYS,
    np,
):
    def calculateGainRate(df, qPCRtarget):
        # Make pairs of (Up day, Down day)
        up_down_pair = list(zip(UP_DAYS[1:], DOWN_DAYS[:]))
        # print(up_down_pair)

        # Process each plasmid
        for plasmid in PlasmidList:
            if plasmid != 'p15A':
                # Process each pair
                for pair in up_down_pair:
                    # Selection Day result
                    numerator = df[(df['Day'] == pair[0]) & (df['PCN'] == plasmid)][qPCRtarget]
                    # The last day of no selection After the selection day
                    denominator = df[(df['Day'] == pair[1]) & (df['PCN'] == plasmid)][qPCRtarget]
                    # Ensure the indexes align
                    denominator = denominator.reset_index(drop=True)
                    numerator = numerator.reset_index(drop=True)
                    gain = numerator/denominator       
                    gain = gain.dropna()

                    indices = df[(df['Day'] == pair[0]) & (df['PCN'] == plasmid)].index

                    # Assign the calculated gain back
                    df.loc[indices, 'gain'] = gain.values

                    # Take mean of all days  
                    df.loc[indices, 'mean of gain'] = np.mean(df.loc[indices, 'gain'])

                    # Calculate the rate & stats
                    generation = np.log2(TreatmentDilution^SelectionTreatmentDays)    
                    gain_rate = np.log(gain)/generation
                    gain_rate_mean = np.mean(gain_rate)
                    gain_rate_std = np.std(gain_rate)

                    # Assign calculated gain rate & stats values to the df
                    df.loc[indices, 'rate_of_gain_log('+qPCRtarget+')'] = gain_rate.values
                    df.loc[indices, 'mean_of_gain_rate(log('+qPCRtarget+'))'] = gain_rate_mean
                    df.loc[indices, 'std_of_gain_rate(log('+qPCRtarget+'))'] = gain_rate_std
                    df.loc[indices, 'ci_lower('+qPCRtarget+')'] = gain_rate_mean - gain_rate_std
                    df.loc[indices, 'ci_upper('+qPCRtarget+')'] = gain_rate_mean + gain_rate_std

        df = df.fillna(np.nan)

        return df
    return calculateGainRate,


@app.cell
def __(calculateGainRate, pd):
    def makeGainPlotDFs(currentCompleteDF):
        currentCompleteDF_gainCopy = currentCompleteDF.copy(deep=True)
        # for the transposon result
        currentCompleteDF_gain_transposon_df = calculateGainRate(currentCompleteDF_gainCopy, 'transposon copy')
        # for the plasmid result
        currentCompleteDF_gain_plasmid_df = calculateGainRate(currentCompleteDF_gainCopy, 'plasmid copy')
        # for the T/P result
        currentCompleteDF_gain_ratio_df = calculateGainRate(currentCompleteDF_gainCopy, 'T/P')

        # To postprocess for ploting
        currentCompleteDF_gain_transposon_df_cleaned_scatter = currentCompleteDF_gain_transposon_df.groupby(['PCN','ExpRep'])[['rate_of_gain_log(transposon copy)']].mean().reset_index().dropna()

        currentCompleteDF_gain_plasmid_df_cleaned_scatter = currentCompleteDF_gain_plasmid_df.groupby(['PCN','ExpRep'])[['rate_of_gain_log(plasmid copy)']].mean().reset_index().dropna()

        currentCompleteDF_gain_ratio_df_cleaned_scatter = currentCompleteDF_gain_ratio_df.groupby(['PCN','ExpRep'])[['rate_of_gain_log(T/P)']].mean().reset_index().dropna()

        # Identify overlapping column names
        overlapping_columns = set(currentCompleteDF_gain_transposon_df_cleaned_scatter.columns) & set(currentCompleteDF_gain_plasmid_df_cleaned_scatter.columns) & set(currentCompleteDF_gain_ratio_df_cleaned_scatter.columns)

        # Concatenate DataFrames horizontally, keeping only one occurrence of overlapping columns
        result_df = pd.concat([currentCompleteDF_gain_transposon_df_cleaned_scatter, currentCompleteDF_gain_plasmid_df_cleaned_scatter.loc[:, ~currentCompleteDF_gain_plasmid_df_cleaned_scatter.columns.isin(overlapping_columns)]], axis=1)

        result_df = pd.concat([result_df, currentCompleteDF_gain_ratio_df_cleaned_scatter.loc[:, ~currentCompleteDF_gain_ratio_df_cleaned_scatter.columns.isin(overlapping_columns)]], axis=1)

        return result_df
    return makeGainPlotDFs,


@app.cell
def __(
    makeGainPlotDFs,
    pd,
    tet10_kan_DF_Use_grouped,
    tet10_noKan_DF_Use_grouped,
):
    scatter_gain_kan = makeGainPlotDFs(tet10_kan_DF_Use_grouped)
    scatter_gain_noKan = makeGainPlotDFs(tet10_noKan_DF_Use_grouped)

    scatter_gain_kan['Kan'] = 1
    scatter_gain_noKan['Kan'] = 0

    combined_scatter_gain_DF = pd.concat([scatter_gain_kan, scatter_gain_noKan])

    combined_scatter_gain_DF['PCN'] = pd.Series(combined_scatter_gain_DF['PCN'], dtype="category")
    combined_scatter_gain_DF["PCN"] = combined_scatter_gain_DF["PCN"].cat.reorder_categories(
                    [
                    'pBR322',
                    'CloDF',
                    'pUC',
                    ],
                ordered=True,
                    )
    return combined_scatter_gain_DF, scatter_gain_kan, scatter_gain_noKan


@app.cell
def __(averagePCN, combined_scatter_DF, combined_scatter_gain_DF):
    combined_scatter_gain_DF.loc[combined_scatter_DF['PCN'] == 'pBR322', 'averagePCN'] = averagePCN[averagePCN['PCN'] == 'pBR322']['plasmid copy'].values[0]
    combined_scatter_gain_DF.loc[combined_scatter_DF['PCN'] == 'CloDF', 'averagePCN'] = averagePCN[averagePCN['PCN'] == 'CloDF']['plasmid copy'].values[0]
    combined_scatter_gain_DF.loc[combined_scatter_DF['PCN'] == 'pUC', 'averagePCN'] = averagePCN[averagePCN['PCN'] == 'pUC']['plasmid copy'].values[0]
    return


@app.cell
def __(alt, colorScheme, combined_scatter_gain_DF):
    # Gain, TCN

    alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title=None),
        alt.Y("rate_of_gain_log(transposon copy)", scale=alt.Scale(domain=[-0.01, 0.24]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, tickCount = 3, ticks=False, grid=False), title=None),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='TCN gain rate',
        width=150,
        height=150
    )
    return


@app.cell
def __():
    return


@app.cell
def __(alt, colorScheme, combined_scatter_gain_DF):
    # Gain, PCN

    alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title=None),
        alt.Y("rate_of_gain_log(plasmid copy)", scale=alt.Scale(domain=[-0.01, 0.24]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, tickCount = 3, ticks=False, grid=False), title=None),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='PCN gain rate',
        width=150,
        height=150
    ).interactive()
    return


@app.cell
def __(combined_scatter_gain_DF):
    combined_scatter_gain_DF
    return


@app.cell
def __(alt, colorScheme, combined_scatter_gain_DF):
    # T/P, gain

    alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title=None),
        alt.Y("rate_of_gain_log(T/P)", scale=alt.Scale(type='log'), axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title=None),
        #domain=[-0.03, 0.33], 
        color=alt.Color('PCN', legend=None).scale(scheme=colorScheme,reverse=True),
       shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='T/P gain rate',
        width=150,
        height=150
    )
    return


@app.cell
def __(alt, colorScheme, combined_scatter_gain_DF):
    TCN_gain = alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_gain_log(transposon copy)", scale=alt.Scale(domain=[0, 0.3]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='gain rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='TCN',
        width=200,
        height=200
    )

    PCN_gain = alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_gain_log(plasmid copy)", scale=alt.Scale(domain=[0, 0.3]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='gain rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='PCN',
        width=200,
        height=200
    )

    Ratio_gain = alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_gain_log(T/P)", scale=alt.Scale(domain=[0, 0.3]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='gain rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
       shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='T/P',
        width=200,
        height=200
    )

    Gain_combined_chart = TCN_gain | PCN_gain | Ratio_gain

    Gain_combined_chart.interactive()
    return Gain_combined_chart, PCN_gain, Ratio_gain, TCN_gain


@app.cell
def __(alt, colorScheme, combined_scatter_gain_DF):
    TCN_gain_PCNasX = alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("averagePCN", axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_gain_log(transposon copy)", scale=alt.Scale(domain=[0, 0.3]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='gain rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='TCN',
        width=200,
        height=200
    )

    PCN_gain_PCNasX = alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("averagePCN", axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_gain_log(plasmid copy)", scale=alt.Scale(domain=[0, 0.3]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='gain rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='PCN',
        width=200,
        height=200
    )

    Ratio_gain_PCNasX = alt.Chart(combined_scatter_gain_DF).mark_point(size=70).encode(
        alt.X("averagePCN", axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_gain_log(T/P)", scale=alt.Scale(domain=[0, 0.3]), axis=alt.Axis(labelFontSize=30, titleFontSize=15, grid=False), title='gain rate'),
        color=alt.Color('PCN').scale(scheme=colorScheme,reverse=True),
        shape=alt.Shape('Kan:O', legend=None, scale=alt.Scale(range=['arrow', 'triangle']))
    ).properties(
        title='T/P',
        width=200,
        height=200
    )

    Gain_combined_PCNasX_chart = TCN_gain_PCNasX | PCN_gain_PCNasX | Ratio_gain_PCNasX

    Gain_combined_PCNasX_chart.interactive()
    return (
        Gain_combined_PCNasX_chart,
        PCN_gain_PCNasX,
        Ratio_gain_PCNasX,
        TCN_gain_PCNasX,
    )


@app.cell
def __(combined_scatter_DF, combined_scatter_gain_DF):
    combined_scatter_gain_DF['Jumping'] = 1
    combined_scatter_DF['Jumping'] = 1
    return


@app.cell
def __(combined_scatter_DF, combined_scatter_gain_DF):
    combined_scatter_gain_DF.to_csv('withJumping_gain_qPCR.csv', index=False)
    combined_scatter_DF.to_csv('withJumping_loss_qPCR.csv', index=False)
    return


@app.cell
def __(combined_data):
    # Export raw data too
    combined_data.to_csv('withJumping_rawData.csv', index=False)
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
