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

    from functools import reduce

    # plotting packages
    import altair as alt

    # Import color utilities from matplotlib
    from statannotations.Annotator import Annotator
    from matplotlib.colors import to_rgba
    from matplotlib.legend_handler import HandlerTuple
    import matplotlib.pyplot as plt
    import seaborn as sns
    return (
        Annotator,
        HandlerTuple,
        alt,
        mo,
        np,
        pd,
        pl,
        plt,
        reduce,
        sns,
        to_rgba,
    )


@app.cell
def __(pd):
    # Import GFP files
    noJumping_gain_df = pd.read_csv('noJumping_gain.csv')
    noJumping_loss_df = pd.read_csv('noJumping_loss.csv')

    Jumping_gain_df = pd.read_csv('withJumping_gain.csv')
    Jumping_loss_df = pd.read_csv('withJumping_loss.csv')

    # Import qPCR files
    noJumping_gain_qPCR_df = pd.read_csv('noJumping_gain_qPCR.csv')
    noJumping_loss_qPCR_df = pd.read_csv('noJumping_loss_qPCR.csv')

    Jumping_gain_qPCR_df = pd.read_csv('withJumping_gain_qPCR.csv')
    Jumping_loss_qPCR_df = pd.read_csv('withJumping_loss_qPCR.csv')
    return (
        Jumping_gain_df,
        Jumping_gain_qPCR_df,
        Jumping_loss_df,
        Jumping_loss_qPCR_df,
        noJumping_gain_df,
        noJumping_gain_qPCR_df,
        noJumping_loss_df,
        noJumping_loss_qPCR_df,
    )


@app.cell
def __(pd):
    # A function to sort dataframe based on plasmid copy number: pBR < CloDF < pUC
    def sortName(df):
        if 'Plasmid' in df.columns.values:       
            df['Plasmid'] = pd.Series(
                df['Plasmid'], dtype="category"
            )  

            if 'Chromosome' in df['Plasmid'].values:
                df['Plasmid'] = df[
                    'Plasmid'
                ].cat.reorder_categories(
                    [
                        "Chromosome",
                        "pBR322",
                        "CloDF13",
                        "pUC",
                    ],
                    ordered=True,       
                )
                # Change Chromosome to Chr
                df['Plasmid'] = df['Plasmid'].replace({
                    'Chromosome': 'Chr',
                })
            else:
                df['Plasmid'] = df[
                    'Plasmid'
                ].cat.reorder_categories(
                    [
                        "pBR322",
                        "CloDF13",
                        "pUC",
                    ],
                    ordered=True,       
            )

            # Change pBR322 to pBR & CloDF13 to CloDF here for consistency
            df['Plasmid'] = df['Plasmid'].replace({
                'pBR322': 'pBR',
                'CloDF13': 'CloDF'
            })

        else:
            df['PCN'] = pd.Series(
                df['PCN'], dtype="category"
            )  
            if 'Chromosome' in df['PCN'].values:
                df['PCN'] = df[
                    'PCN'
                ].cat.reorder_categories(
                    [
                        "Chromosome",
                        "pBR322",
                        "CloDF",
                        "pUC",
                    ],
                    ordered=True,
                )

                # Change Chromosome to Chr
                df['PCN'] = df['PCN'].replace({
                    'Chromosome': 'Chr',
                })

            else:
                df['PCN'] = df[
                    'PCN'
                ].cat.reorder_categories(
                    [
                        "pBR322",
                        "CloDF",
                        "pUC",
                    ],
                    ordered=True,
                )

            # Change pBR322 to pBR
            df['PCN'] = df['PCN'].replace({
                'pBR322': 'pBR',
            })

        return df
    return sortName,


@app.cell
def __(
    Jumping_gain_df,
    Jumping_gain_qPCR_df,
    Jumping_loss_df,
    Jumping_loss_qPCR_df,
    noJumping_gain_df,
    noJumping_gain_qPCR_df,
    noJumping_loss_df,
    noJumping_loss_qPCR_df,
    sortName,
):
    ## Rename & add some column names for consistency here:
    # GFP data
    noJumping_gain_df.rename(columns={'rate_of_gain_log(GFP/OD)': 'rate_log(GFP/OD)'}, inplace=True)
    noJumping_loss_df.rename(columns={'rate_of_loss_log(GFP/OD)': 'rate_log(GFP/OD)'}, inplace=True)
    Jumping_gain_df.rename(columns={'rate_of_gain_log(GFP/OD)': 'rate_log(GFP/OD)'}, inplace=True)
    Jumping_loss_df.rename(columns={'rate_of_loss_log(GFP/OD)': 'rate_log(GFP/OD)'}, inplace=True)

    noJumping_gain_df.rename(columns={'Jumping': 'Transposase'}, inplace=True)
    noJumping_loss_df.rename(columns={'Jumping': 'Transposase'}, inplace=True)
    Jumping_gain_df.rename(columns={'Jumping': 'Transposase'}, inplace=True)
    Jumping_loss_df.rename(columns={'Jumping': 'Transposase'}, inplace=True)

    noJumping_gain_df['gain'] = 1
    noJumping_loss_df['gain'] = 0
    Jumping_gain_df['gain'] = 1
    Jumping_loss_df['gain'] = 0

    noJumping_gain_df_Use = noJumping_gain_df[['Plasmid', 'Kan', 'rate_log(GFP/OD)', 'Transposase', 'gain']]
    noJumping_loss_df_Use = noJumping_loss_df[['Plasmid', 'Kan', 'rate_log(GFP/OD)', 'Transposase', 'gain']]
    Jumping_gain_df_Use = Jumping_gain_df[['Plasmid', 'Kan', 'rate_log(GFP/OD)', 'Transposase', 'gain']]
    Jumping_loss_df_Use = Jumping_loss_df[['Plasmid', 'Kan', 'rate_log(GFP/OD)', 'Transposase', 'gain']]

    noJumping_gain_df_Use = sortName(noJumping_gain_df_Use)
    noJumping_loss_df_Use = sortName(noJumping_loss_df_Use)
    Jumping_gain_df_Use = sortName(Jumping_gain_df_Use)
    Jumping_loss_df_Use = sortName(Jumping_loss_df_Use)

    # qPCR data
    noJumping_gain_qPCR_df.rename(columns={'rate_of_gain_log(transposon copy)': 'rate_log(transposon copy)'}, inplace=True)
    noJumping_loss_qPCR_df.rename(columns={'rate_of_loss_log(transposon copy)': 'rate_log(transposon copy)'}, inplace=True)
    Jumping_gain_qPCR_df.rename(columns={'rate_of_gain_log(transposon copy)': 'rate_log(transposon copy)'}, inplace=True)
    Jumping_loss_qPCR_df.rename(columns={'rate_of_loss_log(transposon copy)': 'rate_log(transposon copy)'}, inplace=True)

    noJumping_gain_qPCR_df.rename(columns={'rate_of_gain_log(plasmid copy)': 'rate_log(plasmid copy)'}, inplace=True)
    noJumping_loss_qPCR_df.rename(columns={'rate_of_loss_log(plasmid copy)': 'rate_log(plasmid copy)'}, inplace=True)
    Jumping_gain_qPCR_df.rename(columns={'rate_of_gain_log(plasmid copy)': 'rate_log(plasmid copy)'}, inplace=True)
    Jumping_loss_qPCR_df.rename(columns={'rate_of_loss_log(plasmid copy)': 'rate_log(plasmid copy)'}, inplace=True)

    noJumping_gain_qPCR_df.rename(columns={'rate_of_gain_log(T/P)': 'rate_log(T/P)'}, inplace=True)
    noJumping_loss_qPCR_df.rename(columns={'rate_of_loss_log(T/P)': 'rate_log(T/P)'}, inplace=True)
    Jumping_gain_qPCR_df.rename(columns={'rate_of_gain_log(T/P)': 'rate_log(T/P)'}, inplace=True)
    Jumping_loss_qPCR_df.rename(columns={'rate_of_loss_log(T/P)': 'rate_log(T/P)'}, inplace=True)

    noJumping_gain_qPCR_df.rename(columns={'Jumping': 'Transposase'}, inplace=True)
    noJumping_loss_qPCR_df.rename(columns={'Jumping': 'Transposase'}, inplace=True)
    Jumping_gain_qPCR_df.rename(columns={'Jumping': 'Transposase'}, inplace=True)
    Jumping_loss_qPCR_df.rename(columns={'Jumping': 'Transposase'}, inplace=True)

    noJumping_gain_qPCR_df['gain'] = 1
    noJumping_loss_qPCR_df['gain'] = 0
    Jumping_gain_qPCR_df['gain'] = 1
    Jumping_loss_qPCR_df['gain'] = 0

    noJumping_gain_qPCR_df_Use = noJumping_gain_qPCR_df[['PCN', 'Kan', 'rate_log(transposon copy)', 'rate_log(plasmid copy)', 'rate_log(T/P)', 'Transposase', 'gain']]
    noJumping_loss_qPCR_df_Use = noJumping_loss_qPCR_df[['PCN', 'Kan', 'rate_log(transposon copy)', 'rate_log(plasmid copy)', 'rate_log(T/P)', 'Transposase', 'gain']]
    Jumping_gain_qPCR_df_Use = Jumping_gain_qPCR_df[['PCN', 'Kan', 'rate_log(transposon copy)', 'rate_log(plasmid copy)', 'rate_log(T/P)', 'Transposase', 'gain']]
    Jumping_loss_qPCR_df_Use = Jumping_loss_qPCR_df[['PCN', 'Kan', 'rate_log(transposon copy)', 'rate_log(plasmid copy)', 'rate_log(T/P)', 'Transposase', 'gain']]

    noJumping_gain_qPCR_df_Use = sortName(noJumping_gain_qPCR_df_Use)
    noJumping_loss_qPCR_df_Use = sortName(noJumping_loss_qPCR_df_Use)
    Jumping_gain_qPCR_df_Use = sortName(Jumping_gain_qPCR_df_Use)
    Jumping_loss_qPCR_df_Use = sortName(Jumping_loss_qPCR_df_Use)
    return (
        Jumping_gain_df_Use,
        Jumping_gain_qPCR_df_Use,
        Jumping_loss_df_Use,
        Jumping_loss_qPCR_df_Use,
        noJumping_gain_df_Use,
        noJumping_gain_qPCR_df_Use,
        noJumping_loss_df_Use,
        noJumping_loss_qPCR_df_Use,
    )


@app.cell
def __(
    Jumping_gain_df_Use,
    Jumping_gain_qPCR_df_Use,
    Jumping_loss_df_Use,
    Jumping_loss_qPCR_df_Use,
    noJumping_gain_df_Use,
    noJumping_gain_qPCR_df_Use,
    noJumping_loss_df_Use,
    noJumping_loss_qPCR_df_Use,
    pd,
):
    combined_df = pd.concat([noJumping_gain_df_Use, noJumping_loss_df_Use, Jumping_gain_df_Use, Jumping_loss_df_Use], ignore_index=True)

    # Concatenate the existing DataFrame with the new rows
    # combined_df = pd.concat([combined_df, noResponse_pUC], ignore_index=True)
    # combined_df = sortName(combined_df)

    # qPCR data processing
    combined_qPCR_df = pd.concat([noJumping_gain_qPCR_df_Use, noJumping_loss_qPCR_df_Use, Jumping_gain_qPCR_df_Use, Jumping_loss_qPCR_df_Use], ignore_index=True)

    # Concatenate the existing DataFrame with the new rows
    # combined_qPCR_df = pd.concat([combined_qPCR_df, noResponse_pUC_qPCR], ignore_index=True)
    combined_qPCR_df['rate_log(transposon copy)'] = combined_qPCR_df['rate_log(transposon copy)'].fillna(0)
    # combined_qPCR_df = sortName(combined_qPCR_df)
    return combined_df, combined_qPCR_df


@app.cell
def __(combined_df, pd):
    # Organize data based on plasmid copy number
    combined_df["Plasmid"] = pd.Series(
            combined_df["Plasmid"], dtype="category"
        )

    combined_df["Plasmid"] = combined_df[
        "Plasmid"
    ].cat.reorder_categories(
        [
            "Chr",
            "pBR",
            "CloDF",
            "pUC",    
        ],
        ordered=True,
    )
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"## qPCR ~ GFP/OD correlation")
    return


@app.cell
def __(pd):
    # GFP files
    noJumping_GFP_df = pd.read_csv('noJumping_rawData_GFP.csv')
    Jumping_GFP_df = pd.read_csv('withJumping_rawData_GFP.csv')

    # qPCR files
    noJumping_qPCR_df = pd.read_csv('noJumping_rawData.csv')
    Jumping_qPCR_df = pd.read_csv('withJumping_rawData.csv')

    noJumping_qPCR_df = noJumping_qPCR_df.rename(columns={'PCN': 'Plasmid'})
    Jumping_qPCR_df = Jumping_qPCR_df.rename(columns={'PCN': 'Plasmid'})

    noJumping_qPCR_df = noJumping_qPCR_df.rename(columns={'ExpRep': 'Replicate'})
    Jumping_qPCR_df = Jumping_qPCR_df.rename(columns={'ExpRep': 'Replicate'})
    return (
        Jumping_GFP_df,
        Jumping_qPCR_df,
        noJumping_GFP_df,
        noJumping_qPCR_df,
    )


@app.cell
def __(Jumping_qPCR_df, noJumping_qPCR_df):
    # Replace multiple values for plotting consistency
    noJumping_qPCR_df['Plasmid'] = noJumping_qPCR_df['Plasmid'].replace({'CloDF': 'CloDF13'})
    Jumping_qPCR_df['Plasmid'] = Jumping_qPCR_df['Plasmid'].replace({'CloDF': 'CloDF13'})

    noJumping_qPCR_df['Kan'] = noJumping_qPCR_df['Kan'].replace({'Kan': 1})
    Jumping_qPCR_df['Kan'] = Jumping_qPCR_df['Kan'].replace({'Kan': 1})
    noJumping_qPCR_df['Kan'] = noJumping_qPCR_df['Kan'].replace({'noKan': 0})
    Jumping_qPCR_df['Kan'] = Jumping_qPCR_df['Kan'].replace({'noKan': 0})
    return


@app.cell
def __(noJumping_GFP_df, noJumping_qPCR_df, pd):
    # Merge the two DataFrames on matching columns
    merged_df = pd.merge(noJumping_qPCR_df, noJumping_GFP_df, on=['Plasmid', 'Day', 'Replicate', 'Kan'])
    return merged_df,


@app.cell
def __(alt, merged_df):
    # Scatter plot of GFP/OD vs transposon copy
    alt.Chart(merged_df).mark_point(size=100, opacity=0.2).encode(
        alt.X("log10(GFP/OD):Q", scale=alt.Scale(domain=[2.8, 5]), axis=alt.Axis(labelFontSize=30, tickCount=2, titleFontSize=15, grid=False), title= None),
        alt.Y("log(transposon copy):Q", scale=alt.Scale(domain=[-5, 5]), axis=alt.Axis(labelFontSize=30, tickCount=2, ticks=False, grid=False, title=None)),
        # 
    ).properties(
        title='',
        width=150,
        height=300
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).interactive()
    return


@app.cell
def __(alt, merged_df):
    # Scatter plot of GFP/OD vs transposon copy
    alt.Chart(merged_df).mark_point(size=100, opacity=0.2).encode(
        alt.X("log10(GFP/OD):Q", scale=alt.Scale(domain=[2.8, 5]), axis=alt.Axis(labelFontSize=30, values = [3, 5], titleFontSize=15, grid=False), title= None),
        alt.Y("log(transposon copy):Q", scale=alt.Scale(domain=[0, 5]), axis=alt.Axis(labelFontSize=30, tickCount = 1, ticks=False, grid=False, title=None)),
        # 
    ).properties(
        title='',
        width=150,
        height=150
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).interactive()
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"# Fig 3 & 4 & all their supplementary figures")
    return


@app.cell
def __(Annotator, HandlerTuple, plt, sns):
    ### Function to barplot qPCR data w/ stats test
    def plotStatQPCR(df, hueVar='Kan', ylim=0.5, yValToPlot="rate_log(transposon copy)"):
        # Create figure with appropriate size
        plt.figure(figsize=(5, 3))

        palette = ["#FAE5D4", "#EA9953", "#DDA8BB", "#BC5078", "#DCC6E7", "#7316A2", ] # "#a50026", 

        # Define the pairs to compare with statistical tests
        # For example, comparing the two Kan groups within each PCN category
        pairs = [
            (("pBR", 1), ("pBR", 0)),  # Replace 0 and 1 with your actual Kan values
            (("CloDF", 1), ("CloDF", 0)),
            (("pUC", 1), ("pUC", 0))
        ]

        # Check if transposons on chromosome only strain is in the data
        if 'Chr' in df.PCN.values:   
            # Create the grouped barplot
            ax = sns.barplot(
                data=df,
                x="PCN",
                y=yValToPlot,
                order=["Chr", "pBR", "CloDF", 'pUC'],
                hue=hueVar,
                # hue_order=[1, 0],
                errorbar="sd",
                capsize=0.1,
                palette=palette,  # Use our custom palette
                # edgecolor='black',  
                # linewidth=1.0,  
                legend=False 
            )
            
            # Assign values to proper legend
            if ax.get_legend() is not None:
                for bars, colors in zip(ax.containers, (["#D5D5D5"]+palette[0::2], palette[1::2])):
                     for bar, color in zip(bars, colors):
                         bar.set_facecolor(color)
                         ax.legend(handles=[tuple(bar_group) for bar_group in ax.containers], labels=[bar_group.get_label() for bar_group in ax.containers], title=ax.legend_.get_title().get_text(), handlelength=4, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.1)})
            else:
                for bars, colors in zip(ax.containers, (["#D5D5D5"]+palette[0::2], palette[1::2])):
                     for bar, color in zip(bars, colors):
                         bar.set_facecolor(color)
                         
            # Aesthetic option: no boarder due to crowdedness
            sns.despine(left=True, bottom=True, right=True, top=True)

            # Create the annotator object
            annotator = Annotator(
                ax, 
                pairs, 
                data=df,
                x="PCN", 
                y=yValToPlot, 
                hue=hueVar,
                order=["Chr", "pBR", "CloDF", 'pUC']
            )

        else:        
            # Create the grouped barplot
            ax = sns.barplot(
                data=df,
                x="PCN",
                y=yValToPlot,
                order=["pBR", "CloDF", 'pUC'],
                hue=hueVar,
                # hue_order=[1, 0],
                errorbar="sd",
                capsize=0.1,
                palette=palette,  # Use our custom palette
                # edgecolor='black',  
                # linewidth=1.0,  
                legend=False 
            )

            # Pair bars and colors
            if ax.get_legend() is not None:
                for bars, colors in zip(ax.containers, (palette[0::2], palette[1::2])):
                     for bar, color in zip(bars, colors):
                         bar.set_facecolor(color)
                         ax.legend(handles=[tuple(bar_group) for bar_group in ax.containers], labels=[bar_group.get_label() for bar_group in ax.containers], title=ax.legend_.get_title().get_text(), handlelength=4, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.1)})
            else:
                for bars, colors in zip(ax.containers, (palette[0::2], palette[1::2])):
                     for bar, color in zip(bars, colors):
                         bar.set_facecolor(color)

            # Aesthetic option: no boarder due to crowdedness
            sns.despine(left=True, bottom=True, right=True, top=True)

            # Create the annotator object
            annotator = Annotator(
                ax, 
                pairs, 
                data=df,
                x="PCN", 
                y=yValToPlot, 
                hue=hueVar,
                order=["pBR", "CloDF", 'pUC']
            )

        # Perform statistical tests and add annotations
        # You can choose different statistical tests like "t-test_ind", "Mann-Whitney", etc.
        annotator.configure(test='Mann-Whitney', text_format='star', loc='outside', fontsize=20)
        annotator.apply_and_annotate()

        # Rest of your code for customizing the plot remains the same
        plt.ylabel(r'$\mathrm{gen}^{-1}$', fontsize=30)
        plt.xlabel("")
        plt.ylim(-0.01, ylim)
        # ax.set_xticklabels([])
        ax.set_yticks([0, ylim])
        ax.set_yticklabels(['0', ylim])
        ax.tick_params(axis='y', labelsize=25,  left=False)
        ax.tick_params(axis='x', labelsize=25)

        ## Uncomment below if want to see legend
        # plt.legend(
        #     title='Kan',
        #     loc='upper left',
        #     fontsize=14,
        #     title_fontsize=16,
        #     markerscale=2,
        #     frameon=True,
        #     handlelength=3,
        #     handleheight=1.5
        # )

        plt.tight_layout()
        plt.show()
    return plotStatQPCR,


@app.cell
def __(Annotator, HandlerTuple, plt, sns):
    ### Function to barplot GFP data w/ stats test
    def plotStatGFP(df, hueVar='Kan', ylim=0.5):
        # Create figure with appropriate size
        plt.figure(figsize=(5, 3))

        palette = ["#FAE5D4", "#EA9953", "#DDA8BB", "#BC5078", "#DCC6E7", "#7316A2", ] # "#a50026", 

        pairs  = [
        (("pBR", 1), ("pBR", 0)),  # 0 and 1 for Kan/Transposase Values
        (("CloDF", 1), ("CloDF", 0)),
        (("pUC", 1), ("pUC", 0))
        ]

        # Check if transposons on chromosome only strain is in the data
        if 'Chr' in df.Plasmid.values:    
            # Create the grouped barplot
            ax = sns.barplot(
                data=df,
                x="Plasmid",
                y="rate_log(GFP/OD)",
                order=["Chr", "pBR", "CloDF", 'pUC'],
                hue=hueVar,
                # hue_order=[1, 0],
                errorbar="sd",
                capsize=0.1,
                palette=palette,  # Use our custom palette
                # edgecolor='black',  
                # linewidth=1.0,     
                legend=False 
            )
            
            # Assign values to proper legend
            if ax.get_legend() is not None:
                for bars, colors in zip(ax.containers, (["#D5D5D5"]+palette[0::2], palette[1::2])):
                     for bar, color in zip(bars, colors):
                         bar.set_facecolor(color)
                         ax.legend(handles=[tuple(bar_group) for bar_group in ax.containers], labels=[bar_group.get_label() for bar_group in ax.containers], title=ax.legend_.get_title().get_text(), handlelength=4, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.1)})
            else:
                for bars, colors in zip(ax.containers, (["#D5D5D5"]+palette[0::2], palette[1::2])):
                     for bar, color in zip(bars, colors):
                         bar.set_facecolor(color)

            # Aesthetic option: no boarder due to crowdedness
            sns.despine(left=True, bottom=True, right=True, top=True)

            # Create the annotator object
            annotator = Annotator(
                ax, 
                pairs,
                data=df,
                x="Plasmid", 
                y="rate_log(GFP/OD)", 
                hue=hueVar,
                order=["Chr","pBR", "CloDF", 'pUC']
            )

        else:
            # Create the grouped barplot
            ax = sns.barplot(
                data=df,
                x="Plasmid",
                y="rate_log(GFP/OD)",
                order=["pBR", "CloDF", 'pUC'],
                hue=hueVar,
                # hue_order=[1, 0],
                errorbar="sd",
                capsize=0.1,
                palette=palette,  # Use our custom palette
                # edgecolor='black',  
                # linewidth=1.0,     
                legend=False 
            )
            
            # Assign values to proper legend
            if ax.get_legend() is not None:
                for bars, colors in zip(ax.containers, (palette[0::2], palette[1::2])):
                     for bar, color in zip(bars, colors):
                         bar.set_facecolor(color)
                         ax.legend(handles=[tuple(bar_group) for bar_group in ax.containers], labels=[bar_group.get_label() for bar_group in ax.containers], title=ax.legend_.get_title().get_text(), handlelength=4, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.1)})
            else:
                for bars, colors in zip(ax.containers, (palette[0::2], palette[1::2])):
                     for bar, color in zip(bars, colors):
                         bar.set_facecolor(color)

            # Aesthetic option: no boarder due to crowdedness
            sns.despine(left=True, bottom=True, right=True, top=True)

            # Create the annotator object
            annotator = Annotator(
                ax, 
                pairs,
                data=df,
                x="Plasmid", 
                y="rate_log(GFP/OD)", 
                hue=hueVar,
                order=["pBR", "CloDF", 'pUC']
            )

        # Perform statistical tests and add annotations
        # You can choose different statistical tests like "t-test_ind", "Mann-Whitney", etc.
        annotator.configure(test='Mann-Whitney', text_format='star', loc='outside', fontsize=20)
        annotator.apply_and_annotate()

        # Rest of your code for customizing the plot remains the same
        plt.ylabel(r'$\mathrm{gen}^{-1}$', fontsize=30)
        plt.xlabel("")
        plt.ylim(-0.01, ylim)
        # ax.set_xticklabels([])
        ax.set_yticks([0, ylim])

        ax.set_yticklabels(['0', ylim])
        ax.tick_params(axis='y', labelsize=25,  left=False)
        ax.tick_params(axis='x', labelsize=25)

        ## Uncomment below if want to see legend
        # plt.legend(
        #     title='Kan',
        #     loc='upper left',
        #     fontsize=14,
        #     title_fontsize=16,
        #     markerscale=2,
        #     frameon=True,
        #     handlelength=3,
        #     handleheight=1.5
        # )

        plt.tight_layout()
        plt.show()
    return plotStatGFP,


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"# Fig. 3")
    return


@app.cell
def __(combined_df, plotStatGFP):
    filtered_df_loss = combined_df[(combined_df['gain']==0) & (combined_df['Kan']==1)]
    # filtered_df_loss['rate_log(GFP/OD)'] = filtered_df_loss['rate_log(GFP/OD)'].fillna(0)

    plotStatGFP(filtered_df_loss, 'Transposase', 0.3)
    return filtered_df_loss,


@app.cell
def __(filtered_qPCR_df_loss):
    filtered_qPCR_df_loss
    return


@app.cell
def __(combined_qPCR_df, plotStatQPCR):
    filtered_qPCR_df_loss = combined_qPCR_df[(combined_qPCR_df['gain']==0) & (combined_qPCR_df['Kan']==1)]

    plotStatQPCR(filtered_qPCR_df_loss, 'Transposase', 1.5)
    return filtered_qPCR_df_loss,


@app.cell
def __(filtered_df_gain):
    filtered_df_gain
    return


@app.cell
def __(combined_df, plotStatGFP):
    filtered_df_gain = combined_df[(combined_df['gain']==1) & (combined_df['Kan']==1)]

    plotStatGFP(filtered_df_gain, 'Transposase', 0.2)
    return filtered_df_gain,


@app.cell
def __(filtered_qPCR_df_gain):
    filtered_qPCR_df_gain
    return


@app.cell
def __(combined_qPCR_df, plotStatQPCR):
    filtered_qPCR_df_gain = combined_qPCR_df[(combined_qPCR_df['gain']==1) & (combined_qPCR_df['Kan']==1)]

    plotStatQPCR(filtered_qPCR_df_gain, 'Transposase', 0.3)
    return filtered_qPCR_df_gain,


@app.cell
def __():
    return


@app.cell
def __():
    # Supplement
    return


@app.cell
def __(combined_df, plotStatGFP):
    filtered_df_loss_noT = combined_df[(combined_df['gain']==0) & (combined_df['Transposase']==0)]
    filtered_df_loss_noT = filtered_df_loss_noT.loc[~(filtered_df_loss_noT['Plasmid'] == 'Chr')]

    plotStatGFP(filtered_df_loss_noT, 'Kan', 0.3)
    return filtered_df_loss_noT,


@app.cell
def __(combined_df, plotStatGFP):
    filtered_df_gain_noT = combined_df[(combined_df['gain']==1) & (combined_df['Transposase']==0)]
    filtered_df_gain_noT = filtered_df_gain_noT.loc[~(filtered_df_gain_noT['Plasmid'] == 'Chr')]

    plotStatGFP(filtered_df_gain_noT, 'Kan', 0.2)
    return filtered_df_gain_noT,


@app.cell
def __():
    return


@app.cell
def __(combined_qPCR_df, plotStatQPCR):
    filtered_qPCR_df_loss_noT = combined_qPCR_df[(combined_qPCR_df['gain']==0) & (combined_qPCR_df['Transposase']==0)]
    filtered_qPCR_df_loss_noT = filtered_qPCR_df_loss_noT.loc[~(filtered_qPCR_df_loss_noT['PCN'] == 'Chr')]

    plotStatQPCR(filtered_qPCR_df_loss_noT, 'Kan', 2)
    return filtered_qPCR_df_loss_noT,


@app.cell
def __(combined_qPCR_df, plotStatQPCR):
    filtered_qPCR_df_gain_noT = combined_qPCR_df[(combined_qPCR_df['gain']==1) & (combined_qPCR_df['Transposase']==0)]
    filtered_qPCR_df_gain_noT = filtered_qPCR_df_gain_noT.loc[~(filtered_qPCR_df_gain_noT['PCN'] == 'Chr')]

    plotStatQPCR(filtered_qPCR_df_gain_noT, 'Kan', 0.2)
    return filtered_qPCR_df_gain_noT,


@app.cell
def __():
    return


@app.cell
def __(filtered_qPCR_df_loss_noT, plotStatQPCR):
    plotStatQPCR(filtered_qPCR_df_loss_noT, 'Kan', 1, 'rate_log(plasmid copy)')
    return


@app.cell
def __(filtered_qPCR_df_gain_noT, plotStatQPCR):
    plotStatQPCR(filtered_qPCR_df_gain_noT, 'Kan', 0.2, 'rate_log(plasmid copy)')
    return


@app.cell
def __(filtered_qPCR_df_gain_noT):
    filtered_qPCR_df_gain_noT
    return


@app.cell
def __():
    return


@app.cell
def __(combined_qPCR_df, plotStatQPCR):
    filtered_qPCR_df_loss_ = combined_qPCR_df[(combined_qPCR_df['gain']==0) & (combined_qPCR_df['Transposase']==1)]
    filtered_qPCR_df_loss_ = filtered_qPCR_df_loss_.loc[~(filtered_qPCR_df_loss_['PCN'] == 'Chr')]

    plotStatQPCR(filtered_qPCR_df_loss_, 'Kan', 0.4, 'rate_log(T/P)')
    return filtered_qPCR_df_loss_,


@app.cell
def __(combined_qPCR_df, plotStatQPCR):
    filtered_qPCR_df_gain_ = combined_qPCR_df[(combined_qPCR_df['gain']==1) & (combined_qPCR_df['Transposase']==1)]
    filtered_qPCR_df_gain_ = filtered_qPCR_df_gain_.loc[~(filtered_qPCR_df_gain_['PCN'] == 'Chr')]

    plotStatQPCR(filtered_qPCR_df_gain_, 'Kan', 0.4, 'rate_log(T/P)')
    return filtered_qPCR_df_gain_,


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"# Fig. 4")
    return


@app.cell
def __(combined_df, plotStatGFP):
    # Filter your data as in your original code
    filtered_df_gain_kan = combined_df[(combined_df['gain']==1) & (combined_df['Transposase']==1)]
    plotStatGFP(filtered_df_gain_kan, 'Kan', 0.2)
    return filtered_df_gain_kan,


@app.cell
def __(filtered_qpcr_df_gain_kan, filtered_qpcr_df_loss_kan):
    print(filtered_qpcr_df_gain_kan[(filtered_qpcr_df_gain_kan['PCN'] == 'pBR') & (filtered_qpcr_df_gain_kan['Kan'] == 0)]['rate_log(transposon copy)'].mean()/filtered_qpcr_df_gain_kan[(filtered_qpcr_df_gain_kan['PCN'] == 'pBR') & (filtered_qpcr_df_gain_kan['Kan'] == 1)]['rate_log(transposon copy)'].mean())

    print(filtered_qpcr_df_loss_kan[(filtered_qpcr_df_loss_kan['PCN'] == 'pBR') & (filtered_qpcr_df_loss_kan['Kan'] == 0)]['rate_log(transposon copy)'].mean()/filtered_qpcr_df_loss_kan[(filtered_qpcr_df_loss_kan['PCN'] == 'pBR') & (filtered_qpcr_df_loss_kan['Kan'] == 1)]['rate_log(transposon copy)'].mean())


    print(filtered_qpcr_df_gain_kan[(filtered_qpcr_df_gain_kan['PCN'] == 'CloDF') & (filtered_qpcr_df_gain_kan['Kan'] == 0)]['rate_log(transposon copy)'].mean()/filtered_qpcr_df_gain_kan[(filtered_qpcr_df_gain_kan['PCN'] == 'CloDF') & (filtered_qpcr_df_gain_kan['Kan'] == 1)]['rate_log(transposon copy)'].mean())

    print(filtered_qpcr_df_loss_kan[(filtered_qpcr_df_loss_kan['PCN'] == 'CloDF') & (filtered_qpcr_df_loss_kan['Kan'] == 0)]['rate_log(transposon copy)'].mean()/filtered_qpcr_df_loss_kan[(filtered_qpcr_df_loss_kan['PCN'] == 'CloDF') & (filtered_qpcr_df_loss_kan['Kan'] == 1)]['rate_log(transposon copy)'].mean())


    print(filtered_qpcr_df_gain_kan[(filtered_qpcr_df_gain_kan['PCN'] == 'pUC') & (filtered_qpcr_df_gain_kan['Kan'] == 0)]['rate_log(transposon copy)'].mean()/filtered_qpcr_df_gain_kan[(filtered_qpcr_df_gain_kan['PCN'] == 'pUC') & (filtered_qpcr_df_gain_kan['Kan'] == 1)]['rate_log(transposon copy)'].mean())

    print(filtered_qpcr_df_loss_kan[(filtered_qpcr_df_loss_kan['PCN'] == 'pUC') & (filtered_qpcr_df_loss_kan['Kan'] == 0)]['rate_log(transposon copy)'].mean()/filtered_qpcr_df_loss_kan[(filtered_qpcr_df_loss_kan['PCN'] == 'pUC') & (filtered_qpcr_df_loss_kan['Kan'] == 1)]['rate_log(transposon copy)'].mean())
    return


@app.cell
def __(combined_qPCR_df, plotStatQPCR):
    filtered_qpcr_df_gain_kan = combined_qPCR_df[(combined_qPCR_df['gain']==1) & (combined_qPCR_df['Transposase']==1)]
    plotStatQPCR(filtered_qpcr_df_gain_kan, 'Kan', 0.3)
    return filtered_qpcr_df_gain_kan,


@app.cell
def __(filtered_qpcr_df_gain_kan, plotStatQPCR):
    # PCN
    plotStatQPCR(filtered_qpcr_df_gain_kan, 'Kan', 0.3, 'rate_log(plasmid copy)')
    return


@app.cell
def __():
    return


@app.cell
def __(combined_df, plotStatGFP):
    filtered_df_loss_kan = combined_df[(combined_df['gain']==0) & (combined_df['Transposase']==1)]
    plotStatGFP(filtered_df_loss_kan, 'Kan', 0.3)
    return filtered_df_loss_kan,


@app.cell
def __(combined_qPCR_df, plotStatQPCR):
    # Filter your data as in your original code
    filtered_qpcr_df_loss_kan = combined_qPCR_df[(combined_qPCR_df['gain']==0) & (combined_qPCR_df['Transposase']==1)]
    plotStatQPCR(filtered_qpcr_df_loss_kan, 'Kan', 0.4)
    return filtered_qpcr_df_loss_kan,


@app.cell
def __(filtered_qpcr_df_loss_kan, plotStatQPCR):
    plotStatQPCR(filtered_qpcr_df_loss_kan, 'Kan', 0.4, 'rate_log(plasmid copy)')
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
