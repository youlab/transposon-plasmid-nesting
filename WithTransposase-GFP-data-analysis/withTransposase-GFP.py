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
    import matplotlib.pyplot as plt 
    import seaborn as sns
    return alt, mo, np, pd, pl, plt, reduce, sns


@app.cell
def __():
    ## excel format and the fluoresence gain one chose
    ## Check this for each exp, can vary slightly according to platereader setting
    OD_START_ROW = 28
    GFP_START_ROW = 92

    ## The end of each treatment course
    UP_DAYS = [0, 7, 14, 21]
    DOWN_DAYS = [5, 12, 19]

    TreatmentDilution = 30000
    LossTreatmentDays = 5
    SelectionTreatmentDays = 2
    PlasmidList = ["p15A", "pBR322", "CloDF13", "pUC", "chromosome"]

    colorScheme = 'plasma'
    return (
        DOWN_DAYS,
        GFP_START_ROW,
        LossTreatmentDays,
        OD_START_ROW,
        PlasmidList,
        SelectionTreatmentDays,
        TreatmentDilution,
        UP_DAYS,
        colorScheme,
    )


@app.cell
def __():
    # Assign files to their corresponding dats
    Day0 = '0730tetA_tube_D0.xlsx'
    Day1 = '0731tetA_tube_D1.xlsx'
    Day2 = '0801tetA_tube_D2.xlsx'
    Day3 = '0802tetA_tube_D3.xlsx'
    Day4 = '0803tetA_tube_D4.xlsx'
    Day5 = '0804tetA_tube_D5.xlsx'
    Day6 = '0805tetA_tube_D6.xlsx'
    Day7 = '0806tetA_tube_D7.xlsx'
    Day8 = '0807tetA_tube_D8.xlsx'
    Day9 = '0808tetA_tube_D9.xlsx'
    Day10 = '0809tetA_tube_D10.xlsx'
    Day11 = '0810tetA_tube_D11.xlsx'
    Day12 = '0811tetA_tube_D12.xlsx'
    Day13 = '0812tetA_tube_D13.xlsx'
    Day14 = '0813tetA_tube_D14.xlsx'
    Day15 = '0814tetA_tube_D15.xlsx'
    Day16 = '0815tetA_tube_D16.xlsx'
    Day17 = '0816tetA_tube_D17.xlsx'
    Day18 = '0817tetA_tube_D18.xlsx'
    Day19 = '0818tetA_tube_D19.xlsx'
    Day20 = '0819tetA_tube_D20.xlsx'
    Day21 = '0820tetA_tube_D21.xlsx'

    datafiles = [Day0, Day1, Day2, Day3, Day4, Day5, Day6, Day7, Day8, Day9, Day10, Day11, Day12, Day13, Day14, Day15, Day16, Day17, Day18, Day19, Day20, Day21]
    return (
        Day0,
        Day1,
        Day10,
        Day11,
        Day12,
        Day13,
        Day14,
        Day15,
        Day16,
        Day17,
        Day18,
        Day19,
        Day2,
        Day20,
        Day21,
        Day3,
        Day4,
        Day5,
        Day6,
        Day7,
        Day8,
        Day9,
        datafiles,
    )


@app.cell
def __(pd):
    noJumping_rawData_GFP_df = pd.read_csv('noJumping_rawData_GFP.csv')
    return noJumping_rawData_GFP_df,


@app.cell
def __(pd):
    def prepare_metadata(plate_metadata_path):
        metadata_df = pd.read_csv(plate_metadata_path, index_col='<>').rename_axis('PlateRow').reset_index().melt(id_vars=['PlateRow'], var_name='PlateColumn', value_name='Sample')

        metadata_df['Well'] = metadata_df['PlateRow'].astype(str) + metadata_df['PlateColumn'].astype(str)

        metadata_df = metadata_df.drop(columns=['PlateRow', 'PlateColumn'])

        return metadata_df
    return prepare_metadata,


@app.cell
def __(pd):
    def duplicate_and_modify(df, column, value, new_values):
        # Find rows that match the condition
        mask = df[column] == value

        # Create a list to store the new rows
        new_rows = []

        # Iterate through matching rows
        for idx, row in df[mask].iterrows():
            # # Add the original row
            # new_rows.append(row)

            # Create duplicates with modified values
            for new_value in new_values:
                new_row = row.copy()
                new_row[column] = new_value
                new_rows.append(new_row)

        # Create a new DataFrame with the new rows
        new_df = pd.DataFrame(new_rows)

        # Combine the new DataFrame with the original, excluding the original matching rows
        mask_inverse = df[df[column] != value]
        # print(mask_inverse)
        result = pd.concat([mask_inverse, new_df], ignore_index=True)

        return result
    return duplicate_and_modify,


@app.cell
def __(pd):
    def duplicate_and_modify_rows(df, condition_column, condition_value, modify_column, new_value):
        """
        Duplicates rows where condition_column equals condition_value,
        and changes the value in modify_column for the duplicated rows.

        :param df: The input DataFrame
        :param condition_column: The column to check for the condition
        :param condition_value: The value to match in the condition column
        :param modify_column: The column to modify in the duplicated rows
        :param new_value: The new value to set in the modified column
        :return: The modified DataFrame
        """
        # Create a mask for rows that meet the condition
        mask = df[condition_column] == condition_value

        # Duplicate the rows that meet the condition
        duplicated_rows = df[mask].copy()

        # Modify the specified column in the duplicated rows
        duplicated_rows[modify_column] = new_value

        # Concatenate the original DataFrame with the duplicated and modified rows
        result_df = pd.concat([df, duplicated_rows], ignore_index=True)

        # Sort the DataFrame to keep the duplicated rows next to their originals
        result_df = result_df.sort_index().reset_index(drop=True)

        return result_df
    return duplicate_and_modify_rows,


@app.cell
def __(
    Day0,
    GFP_START_ROW,
    OD_START_ROW,
    duplicate_and_modify,
    np,
    pd,
    prepare_metadata,
):
    def prepare_plate_data(excel_filepath, sheetName, day):
        ## import plate layout metadata
        ## Special treatment for tet10, since it does not have the psc101 data
        if excel_filepath == Day0:
            # print('here')
            metadata_df = prepare_metadata("plate-layout-metadata_Day0.csv")

            ## OD600 data
            my_OD_df = (
                pd.read_excel(
                    excel_filepath,
                    sheet_name=sheetName,
                    skiprows=OD_START_ROW - 1,
                    nrows=4,
                    usecols=list(range(0, 11)),
                    index_col="<>",
                )
                .rename_axis("PlateRow")
                .reset_index()
                .melt(
                    id_vars=["PlateRow"], var_name="PlateColumn", value_name="OD"
                )
            )

            ## GFP data
            my_GFP_df = (
                pd.read_excel(
                    excel_filepath,
                    sheet_name=sheetName,
                    skiprows=88 - 1,
                    nrows=4,
                    usecols=list(range(0, 11)),
                    index_col="<>",
                )
                .rename_axis("PlateRow")
                .reset_index()
                .melt(
                    id_vars=["PlateRow"], var_name="PlateColumn", value_name="GFP"
                )
            )

        else:
            if sheetName == "tet10":
                metadata_df = prepare_metadata("plate-layout-metadata_tet10.csv")
            else:
                metadata_df = prepare_metadata("plate-layout-metadata.csv")

            ## OD600 data
            my_OD_df = (
                pd.read_excel(
                    excel_filepath,
                    sheet_name=sheetName,
                    skiprows=OD_START_ROW - 1,
                    nrows=6,
                    usecols=list(range(0, 11)),
                    index_col="<>",
                )
                .rename_axis("PlateRow")
                .reset_index()
                .melt(
                    id_vars=["PlateRow"], var_name="PlateColumn", value_name="OD"
                )
            )

            ## GFP data
            my_GFP_df = (
                pd.read_excel(
                    excel_filepath,
                    sheet_name=sheetName,
                    skiprows=GFP_START_ROW - 1,
                    nrows=6,
                    usecols=list(range(0, 11)),
                    index_col="<>",
                )
                .rename_axis("PlateRow")
                .reset_index()
                .melt(
                    id_vars=["PlateRow"], var_name="PlateColumn", value_name="GFP"
                )
            )

        ## merge OD and GFP data.
        my_plate_df = pd.merge(my_OD_df, my_GFP_df, on=["PlateRow", "PlateColumn"])

        ## make a 'Well' column
        my_plate_df["Well"] = my_plate_df["PlateRow"].astype(str) + my_plate_df[
            "PlateColumn"
        ].astype(str)
        ## drop the 'PlateRow' and 'PlateColumn' columns
        my_plate_df = my_plate_df.drop(columns=["PlateRow", "PlateColumn"])

        ## merge with metadata based on the 'Well' column
        my_plate_df_with_samples = pd.merge(my_plate_df, metadata_df, on=["Well"])

        ## Calculate the mean GFP & OD for LB only
        LB_ref = my_plate_df_with_samples[
            my_plate_df_with_samples["Sample"] == "LBRef"
        ]
        # print(LB_ref)
        LB_ref_GFP = np.mean(LB_ref["GFP"])
        LB_ref_OD = np.mean(LB_ref["OD"])

        ## remove ref now
        my_plate_df_with_samples = my_plate_df_with_samples[
            ~my_plate_df_with_samples["Sample"].str.contains("LBRef")
        ]

        ## remove Samples containing Water or LB
        my_plate_df_with_samples = my_plate_df_with_samples[
            ~my_plate_df_with_samples["Sample"].str.contains("empty|noUse")
        ]

        ## Splitting the 'Sample' column into two new columns
        my_plate_df_with_samples[["Kan_Plasmid", "Replicate"]] = (
            my_plate_df_with_samples["Sample"].str.split("_", expand=True)
        )

        ## Splitting the 'Sample' column again based on with Kan or not
        my_plate_df_with_samples[["Kan", "Plasmid"]] = (
            my_plate_df_with_samples["Kan_Plasmid"].str.split("-", expand=True)
        )

        my_plate_df_with_samples = my_plate_df_with_samples.drop(['Kan_Plasmid', 'Sample'], axis=1)

        # my_plate_df_with_samples = duplicate_and_modify(my_plate_df_with_samples, 'Kan', 'noKan', ['Kan'])

        my_plate_df_with_samples = duplicate_and_modify(my_plate_df_with_samples, 'Replicate', '0', ['1','2','3','4','5'])
        # my_plate_df_with_samples = duplicate_and_modify_rows(my_plate_df_with_samples, 'Plasmid', 'Chromosome', 'Kan', 'Kan')
        # print(my_plate_df_with_samples)

        ## Sort Plasmid by copy number by encoding the Plasmid Column as Categorical Data
        ## https://pandas.pydata.org/docs/user_guide/categorical.html#reordering
        my_plate_df_with_samples["Plasmid"] = pd.Series(
            my_plate_df_with_samples["Plasmid"], dtype="category"
        )
        # print(my_plate_df_with_samples)

        # now reorder the column based on copy number
        ## Check this
        if sheetName == "tet10":
            my_plate_df_with_samples["Plasmid"] = my_plate_df_with_samples[
                "Plasmid"
            ].cat.reorder_categories(
                [
                    "p15A",
                    "pBR322",
                    "CloDF13",
                    "pUC",
                    "Chromosome",

                ],
                ordered=True,
            )
        else:
            my_plate_df_with_samples["Plasmid"] = my_plate_df_with_samples[
                "Plasmid"
            ].cat.reorder_categories(
                [
                    "pSC101",
                    "p15A",
                    "pBR322",
                    "CloDF13",
                    "pUC",
                    "Chromosome",
                ],
                ordered=True,
            )

        ## add the Day of sample collection
        my_plate_df_with_samples["Day"] = day

        # print('hi', my_plate_df_with_samples)

        ## add the PCN of each plasmid
        my_plate_df_with_samples.loc[
            my_plate_df_with_samples["Plasmid"] == "chromosome", "PCN"
        ] = 0

        # my_plate_df_with_samples.loc[
        #     my_plate_df_with_samples["Plasmid"] == "pSC101", "PCN"
        # ] = 5
        my_plate_df_with_samples.loc[
            my_plate_df_with_samples["Plasmid"] == "p15A", "PCN"
        ] = 20
        my_plate_df_with_samples.loc[
            my_plate_df_with_samples["Plasmid"] == "pBR322", "PCN"
        ] = 60
        my_plate_df_with_samples.loc[
            my_plate_df_with_samples["Plasmid"] == "CloDF13", "PCN"
        ] = 100
        my_plate_df_with_samples.loc[
            my_plate_df_with_samples["Plasmid"] == "pUC", "PCN"
        ] = 200

        ## calculate GFP/OD
        my_plate_df_with_samples["GFP/OD"] = (
            my_plate_df_with_samples["GFP"] - LB_ref_GFP
        ) / (my_plate_df_with_samples["OD"] - LB_ref_OD)

        ## calculate log10 (GFP/OD)
        my_plate_df_with_samples["log10(GFP/OD)"] = np.log10(
            my_plate_df_with_samples["GFP/OD"]
        )

        return my_plate_df_with_samples
    return prepare_plate_data,


@app.cell
def __(datafiles, prepare_plate_data):
    # Creat a list of dfs for each treatment, each list contains the df of one day
    dataframes_tet10 = [prepare_plate_data(df, 'tet10', zero_indexed_day) for zero_indexed_day, df in enumerate(datafiles)]
    return dataframes_tet10,


@app.cell
def __(dataframes_tet10, pd):
    completeDF_tet10 = pd.concat(dataframes_tet10, ignore_index=True)
    completeDF_tet10['tet'] = 10
    completeDF_tet10.loc[completeDF_tet10['Kan'] == 'Kan', 'Kan'] = 1
    completeDF_tet10.loc[completeDF_tet10['Kan'] == 'noKan', 'Kan'] = 0

    completeDF_tet10['Selection'] = completeDF_tet10['Day'].isin([0, 6, 7, 13, 14, 20, 21])
    return completeDF_tet10,


@app.cell
def __(completeDF_tet10):
    completeDF_tet10_noP15A = completeDF_tet10[completeDF_tet10['Plasmid'] != 'p15A']
    completeDF_tet10_noP15A_copy = completeDF_tet10_noP15A.copy(deep = True)
    completeDF_tet10_noP15A = completeDF_tet10_noP15A[completeDF_tet10_noP15A['Plasmid'] != 'Chromosome']

    return completeDF_tet10_noP15A, completeDF_tet10_noP15A_copy


@app.cell
def __(alt, colorScheme, completeDF_tet10_noP15A_copy):
    scatter_plot_log_tet10_alone_line = alt.Chart(completeDF_tet10_noP15A_copy).mark_line(
        point=alt.OverlayMarkDef(size=100, opacity=0.7)).encode(
        x=alt.X('Day:O', axis=alt.Axis(labelFontSize=15, labelAngle=0, ticks=False)),
        y=alt.Y('log10(GFP/OD)', scale=alt.Scale(domain=[2.6, 5]), axis=alt.Axis(labelFontSize=15, ticks=False, grid=False)),
        color=alt.Color('Plasmid').scale(scheme=colorScheme,reverse=True),
        shape='Selection',
        detail='Replicate'
    ).properties(
        width=510,
        height=300
    ).facet(
        column=alt.Column('Kan', header=alt.Header(title=None))
    ).interactive()


    scatter_plot_log_tet10_alone_line.configure_axis(
        # labelFontSize=14,
        titleFontSize=16
    ).configure_header(
        titleFontSize=20,
        labelFontSize=14 
    ).properties(
        title=alt.TitleParams(
            text='tet10',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    )
    return scatter_plot_log_tet10_alone_line,


@app.cell
def __(completeDF_tet10_noP15A, noJumping_rawData_GFP_df, pd):
    noJumping_rawData_GFP_df_chr = noJumping_rawData_GFP_df[noJumping_rawData_GFP_df['Plasmid']=='Chromosome']
    # Reset index if needed
    finalWControl = pd.concat([completeDF_tet10_noP15A, noJumping_rawData_GFP_df_chr], ignore_index=True)
    return finalWControl, noJumping_rawData_GFP_df_chr


@app.cell
def __(finalWControl, pd):
    # Sort DF based on PCN + chr
    finalWControl["Plasmid"] = pd.Series(
            finalWControl["Plasmid"], dtype="category"
        )

    finalWControl["Plasmid"] = finalWControl[
                "Plasmid"
            ].cat.reorder_categories(
                [
                    "pBR322",
                    "CloDF13",
                    "pUC",
                    "Chromosome",
                ],
                ordered=True,
            )
    return


@app.cell
def __(completeDF_tet10_noP15A):
    completeDF_tet10_noP15A
    return


@app.cell
def __(completeDF_tet10_noP15A):
    completeDF_tet10_noP15A_Kan = completeDF_tet10_noP15A[completeDF_tet10_noP15A['Kan'] == 1]
    completeDF_tet10_noP15A_NoKan = completeDF_tet10_noP15A[completeDF_tet10_noP15A['Kan'] == 0]
    return completeDF_tet10_noP15A_Kan, completeDF_tet10_noP15A_NoKan


@app.cell
def __(completeDF_tet10_noP15A_Kan, plt, sns):
    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    palette = ["#EA9953", "#BC5078", "#7316A2", "#D5D5D5"]


    # Create line plot with points
    ax = sns.lineplot(
        data=completeDF_tet10_noP15A_Kan[completeDF_tet10_noP15A_Kan['Plasmid']!='Chromosome'],
        x='Day',
        y='log10(GFP/OD)',
        hue='Plasmid',  
        units='Replicate',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=1,       # Match opacity from Altair
        palette=palette ,
        legend=False     # No legend as in original
    )

    # # Position legend in top right corner with customization
    # plt.legend(
    #     loc='center left',      # Base location 
    #     bbox_to_anchor=(1.01, 0.5),  # Place it outside (x, y) coordinates
    #     frameon=False,          # Remove border
    #     fontsize=18,            # Font size
    #     title_fontsize=18       # Title font size
    # )

    # # Get the handles (visual elements like lines or markers)
    # handles = ax.legend_.legend_handles

    # # Create a new legend with your custom labels
    # plt.legend(handles=handles, 
    #            labels=["pBR", "CloDF", "pUC", "Chr"],
    #            loc='center left',      # Base location 
    #             bbox_to_anchor=(1.01, 0.5),  # Place it outside (x, y) coordinates
    #             frameon=False,          # Remove border
    #             fontsize=18,            # Font size
    #             title_fontsize=18       # Title font size
    # )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(2.5, 5)

    # Configure axis appearance
    ax.tick_params(axis='x', labelsize=25, length=0)
    ax.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([3,4, 5])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('ln(GFP/OD)', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax, palette


@app.cell
def __(completeDF_tet10_noP15A_NoKan, plt, sns):
    # Set up figure size
    plt.figure(figsize=(6.5, 3))

    palette02 = ["#7316A2", "#EA9953", "#BC5078"]

    # Create line plot with points
    ax_noKan = sns.lineplot(
        data=completeDF_tet10_noP15A_NoKan[completeDF_tet10_noP15A_NoKan['Plasmid']!='Chromosome'],
        x='Day',
        y='log10(GFP/OD)',
        hue='Plasmid',  
        units='Replicate',  # Creates separate lines for each ExpRep (detail in Altair)
        estimator=None,  # Don't aggregate data points
        marker='o',      # Add markers on the line
        markersize=12,   # Large marker size
        alpha=1,       # Match opacity from Altair
        palette=palette02 ,
        legend=False     # No legend as in original
    )

    # # Position legend in top right corner with customization
    # plt.legend(
    #     loc='center left',      # Base location 
    #     bbox_to_anchor=(1.01, 0.5),  # Place it outside (x, y) coordinates
    #     frameon=False,          # Remove border
    #     fontsize=18,            # Font size
    #     title_fontsize=18       # Title font size
    # )

    # # Get the handles (visual elements like lines or markers)
    # handles = ax.legend_.legend_handles

    # # Create a new legend with your custom labels
    # plt.legend(handles=handles, 
    #            labels=["pBR", "CloDF", "pUC", "Chr"],
    #            loc='center left',      # Base location 
    #             bbox_to_anchor=(1.01, 0.5),  # Place it outside (x, y) coordinates
    #             frameon=False,          # Remove border
    #             fontsize=18,            # Font size
    #             title_fontsize=18       # Title font size
    # )

    # Set axis limits to match Altair scales
    plt.xlim(-1, 22)
    plt.ylim(2.5, 5)

    # Configure axis appearance
    ax_noKan.tick_params(axis='x', labelsize=25, length=0)
    ax_noKan.tick_params(axis='y', labelsize=25, length=0)
    plt.xticks([0, 5, 10, 15, 20])  # 5 tick marks
    plt.yticks([3,4, 5])  # 3 tick marks
    plt.xlabel('')
    plt.ylabel('ln(GFP/OD)', fontsize=30)
    plt.grid(False)

    # Remove spines to match minimal Altair style
    sns.despine(left=True, bottom=True, right=True, top=True)

    # Use tight layout to optimize spacing
    plt.tight_layout()

    # Display plot
    plt.show()
    return ax_noKan, palette02


@app.cell
def __(mo):
    mo.md(r"# Calculate loss rate")
    return


@app.cell
def __(
    DOWN_DAYS,
    LossTreatmentDays,
    PlasmidList,
    TreatmentDilution,
    UP_DAYS,
    np,
    pd,
):
    def calculateLossRateAvg(df):
        return_df_list = []
        for plasmid in PlasmidList:
            if plasmid in df['Plasmid'].values:
                copy_number = list(set(df[df['Plasmid'] == plasmid]['PCN']))[0]

                # print(plasmid)
                loss01 = \
                (df[(df['Day'] == UP_DAYS[0]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index() / \
                df[(df['Day'] == DOWN_DAYS[0]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index()).reset_index()['GFP/OD']

                KanSeries = df[(df['Day'] == UP_DAYS[1]) & (df['Plasmid'] == plasmid)]['Kan'].reset_index(drop=True)
                # print(KanSeries)

                loss02 = \
                (df[(df['Day'] == UP_DAYS[1]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index() / \
                df[(df['Day'] == DOWN_DAYS[1]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index()).reset_index()['GFP/OD']

                loss03 = \
                (df[(df['Day'] == UP_DAYS[2]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index() / \
                df[(df['Day'] == DOWN_DAYS[2]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index()).reset_index()['GFP/OD']

                loss = (loss01+loss02+loss03)/3

                generation = np.log2(TreatmentDilution^LossTreatmentDays)

                loss_rate = np.log(loss)/generation
                # print(loss_rate)
                loss_rate_mean = np.mean(loss_rate)
                loss_rate_std = np.std(loss_rate)

                loss_rate_DF = pd.DataFrame({'Plasmid': plasmid,\
                                             'PCN': copy_number,
                                             'Kan': KanSeries,
                                             'rate_of_loss_log(GFP/OD)': loss_rate})
                # print(loss_rate_DF)

                return_df_list.append(loss_rate_DF)
                # loss, type(loss), lossRate, type(lossRate)

        # reduce(lambda left, right: pd.merge(left, right, on='ID', how='outer'), dfs)
        return_df = pd.concat(return_df_list)
        return_df['Plasmid'] = pd.Series(return_df['Plasmid'], dtype="category")
        return_df["Plasmid"] = return_df[
                "Plasmid"
            ].cat.reorder_categories(
                [
                    "pBR322",
                    "CloDF13",
                    "pUC",
                ],
                ordered=True,
            )


        return return_df
    return calculateLossRateAvg,


@app.cell
def __(calculateLossRateAvg, completeDF_tet10_noP15A):
    lossRate_tet10_df = calculateLossRateAvg(completeDF_tet10_noP15A)
    lossRate_tet10_df['tet'] = 10
    return lossRate_tet10_df,


@app.cell
def __(alt, colorScheme, lossRate_tet10_df):
    alt.Chart(lossRate_tet10_df).mark_circle(size=60).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_loss_log(GFP/OD)", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title='loss rate'),
        color=alt.Color('Plasmid').scale(scheme=colorScheme,reverse=True)
    ).properties(
        width=200,
        height=200
    ).facet(
        column=alt.Column('Kan', header=alt.Header(title=None))
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).properties(
        title=alt.TitleParams(
            text='loss rate per generation',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    ).interactive()
    return


@app.cell
def __():
    return


@app.cell
def __(mo):
    mo.md(r"# Calculate gain rate")
    return


@app.cell
def __(
    DOWN_DAYS,
    LossTreatmentDays,
    PlasmidList,
    TreatmentDilution,
    UP_DAYS,
    np,
    pd,
):
    def calculateGainRateAvg(df):
        return_df_list = []
        for plasmid in PlasmidList:
            if plasmid in df['Plasmid'].values:
                copy_number = list(set(df[df['Plasmid'] == plasmid]['PCN']))[0]

                # print(plasmid)
                gain01 = \
                (df[(df['Day'] == UP_DAYS[1]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index() / \
                df[(df['Day'] == DOWN_DAYS[0]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index()).reset_index()['GFP/OD']

                KanSeries = df[(df['Day'] == UP_DAYS[1]) & (df['Plasmid'] == plasmid)]['Kan'].reset_index(drop=True)
                # print(KanSeries)

                gain02 = \
                (df[(df['Day'] == UP_DAYS[2]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index() / \
                df[(df['Day'] == DOWN_DAYS[1]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index()).reset_index()['GFP/OD']

                # gain03 = \
                # (df[(df['Day'] == UP_DAYS[3]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index() / \
                # df[(df['Day'] == DOWN_DAYS[2]) & (df['Plasmid'] == plasmid)]['GFP/OD'].reset_index()).reset_index()['GFP/OD']

                gain = (gain01+gain02)/2

                generation = np.log2(TreatmentDilution^LossTreatmentDays)

                gain_rate = np.log(gain)/generation
                # print(gain_rate)
                gain_rate_mean = np.mean(gain_rate)
                gain_rate_std = np.std(gain_rate)

                gain_rate_DF = pd.DataFrame({'Plasmid': plasmid,\
                                             'PCN': copy_number,
                                             'Kan': KanSeries,
                                             'rate_of_gain_log(GFP/OD)': gain_rate})
                # print(gain_rate_DF)

                return_df_list.append(gain_rate_DF)
                # gain, type(gain), gainRate, type(gainRate)

        # reduce(lambda left, right: pd.merge(left, right, on='ID', how='outer'), dfs)
        return_df = pd.concat(return_df_list)
        return_df['Plasmid'] = pd.Series(return_df['Plasmid'], dtype="category")


        return return_df
    return calculateGainRateAvg,


@app.cell
def __(calculateGainRateAvg, completeDF_tet10_noP15A):
    gainRate_tet10_df = calculateGainRateAvg(completeDF_tet10_noP15A)
    gainRate_tet10_df['tet'] = 10
    return gainRate_tet10_df,


@app.cell
def __(alt, colorScheme, gainRate_tet10_df):
    alt.Chart(gainRate_tet10_df).mark_circle(size=60).encode(
        alt.X("PCN", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title='PCN'),
        alt.Y("rate_of_gain_log(GFP/OD)", axis=alt.Axis(labelFontSize=15, titleFontSize=15, grid=False), title='gain rate'),
        color=alt.Color('Plasmid').scale(scheme=colorScheme,reverse=True)
    ).properties(
        width=200,
        height=200
    ).facet(
        column=alt.Column('Kan', header=alt.Header(title=None))
    ).configure_header(
        labelFontSize=20  # Adjust this value to change the subtitle size
    ).properties(
        title=alt.TitleParams(
            text='gain rate per generation',
            fontSize=20,
            fontWeight='bold',
            # color='grey',
            anchor='middle'
        )
    ).interactive()
    return


@app.cell
def __(gainRate_tet10_df, pd):
    gainRate_tet10_df['Plasmid'] = pd.Series(gainRate_tet10_df['Plasmid'], dtype="category")
    gainRate_tet10_df["Plasmid"] = gainRate_tet10_df["Plasmid"].cat.reorder_categories(
                    [
                    'pBR322',
                    'CloDF13',
                    'pUC',
                    ],
                ordered=True,
                    )
    return


@app.cell
def __(gainRate_tet10_df, lossRate_tet10_df):
    lossRate_tet10_df['Jumping'] = 1
    gainRate_tet10_df['Jumping'] = 1
    return


@app.cell
def __(gainRate_tet10_df, lossRate_tet10_df):
    gainRate_tet10_df.to_csv('withJumping_gain.csv', index=False)
    lossRate_tet10_df.to_csv('withJumping_loss.csv', index=False)
    return


@app.cell
def __():
    # completeDF_tet10_noP15A.to_csv('withJumping_rawData_GFP.csv', index=False)
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
