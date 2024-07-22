#!/usr/bin/env python

############################
#### Summary of clustering
#### author: Alan Tracey
#### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
############################


import pandas as pd
import sys
import numpy as np


def calculate_zygosity_confidence(filtered_df):
    # Define target values for each classification
    targets = {
        'Hom WT': 100,
        'Het NHEJ': 50,
        'Hom NHEJ': 0
    }

    # Define how strict the confidence measurement should be
    leniency = {
        'Hom WT': 1,
        'Het NHEJ': 0.5,  # More lenient for Het NHEJ
        'Hom NHEJ': 1
    }

    def get_confidence(row):
        # Assuming columns like 'Reads_WT', 'Reads_Mut', etc., sum these to get total reads
        total_reads = sum([row[col] for col in filtered_df.columns if 'Reads' in col])

        # Calculate the confidence based on classification
        target = targets.get(row['Classification'], None)
        if target is None:
            return None

        difference = abs(row['% Wt'] - target)
        adjusted_difference = difference * leniency.get(row['Classification'], 1)
        confidence = max(0, 1 - (adjusted_difference / 100))

        # Adjust confidence based on total reads
        if total_reads < 3000:
            penalty = (3000 - total_reads) / 3000 * 0.1  # Up to 10% penalty for amplicons with fewer than 3000 reads.  Penalty grows with distance below 3000.
            confidence -= penalty
            confidence = max(0, confidence)  # Ensure confidence doesn't go below 0

        return confidence

    # Apply the confidence calculation to each row in the DataFrame
    filtered_df['Class_Conf'] = filtered_df.apply(get_confidence, axis=1)

    return filtered_df



def parse_edits_csv(df):
    # Calculate total reads per row
    df['Total Reads'] = df[
        ['Wt', 'Template-based', 'Delins', 'Ins_inframe', 'Ins_outframe', 'Dels_inframe', 'Dels_outframe']].sum(
        axis=1)

    # Calculate percentage of wild-type reads
    df['% Wt'] = (df['Wt'] / df['Total Reads'] * 100)

    # Calculate percentage deletions
    df['% Dels'] = (df['Dels_inframe'] + df['Dels_outframe']) / df['Total Reads'] * 100

    # Calculate percentage insertions
    df['% Ins'] = (df['Ins_inframe'] + df['Ins_outframe']) / df['Total Reads'] * 100

    # Calculate percentage delins
    df['% Delins'] = df['Delins'] / df['Total Reads'] * 100

    df['Classification'] = df['% Wt'].apply(classify)


    return df


def classify(wt_percentage):
    if wt_percentage > 80:
        return 'Hom WT'
    elif 40 <= wt_percentage <= 60:
        return 'Het NHEJ'
    elif wt_percentage < 20:
        return 'Hom NHEJ'
    else:
        return 'Ambiguous'



def analyze_clonality(grouped_dels, grouped_ins, edits_df, min_read_threshold):
    """
    Analyzes clonality by examining peak distributions within editing data.

    Parameters:
    - grouped_dels (DataFrame): DataFrame containing grouped deletion data.
    - grouped_ins (DataFrame): DataFrame containing grouped insertion data.
    - edits_df (DataFrame): DataFrame containing edits and associated metrics.
    - min_read_threshold (int): Minimum read count required for valid analysis.

    Returns:
    - dict: Dictionary containing various metrics related to clonality analysis.
    """

    # Check for insufficient data for analysis
    if edits_df.empty or 'Total Reads' not in edits_df.columns or edits_df['Total Reads'].iloc[0] < min_read_threshold:
        return {
            "Classification": "Ambiguous",
            "edition_peak_count": 0,
            "clonality": "Insufficient data for analysis"
        }

    # Combine deletion and insertion data, calculate proportions
    combined_df = pd.concat([grouped_dels, grouped_ins])
    total_counts = edits_df['Total Reads'].iloc[0]
    combined_df['Proportion'] = combined_df['Count'] / total_counts

    # Determine significant peaks
    significant_peaks = combined_df[combined_df['Proportion'] > 0.05]
    peak_proportions = significant_peaks['Proportion'].tolist()

    # Calculate metrics to assess clonality
    max_peak = significant_peaks['Proportion'].max() if not significant_peaks.empty else 0
    wt_perc = edits_df['% Wt'].iloc[0] if not edits_df.empty else 0
    peak_occupancy = sum(significant_peaks['Proportion']) if not significant_peaks.empty else 0

    # Evaluate the distribution and dominance of peaks
    dominant_peak_proportion = max_peak
    sum_of_other_peaks = peak_occupancy - dominant_peak_proportion

    # Clonality categorization logic
    if wt_perc > 85:
        clonality = "Low editing activity"
    elif dominant_peak_proportion > 0.85:
        clonality = "Clonal"
    elif len(significant_peaks) == 1 and max_peak > 0.4 and wt_perc > 0.4:
        clonality = "Clonal"
    elif len(significant_peaks) == 2 and peak_occupancy >= 0.8:
        clonality = "Clonal"
    elif (len(significant_peaks) in [1, 2]) and peak_occupancy > 0.75:
        clonality = "Likely clonal with minor background variants"
    elif len(significant_peaks) > 2 and sum_of_other_peaks > 0.4:
        clonality = "Polyclonal"
    else:
        clonality = "Ambiguous"

    # Re-calculate zygosity confidence for updated clonality categorization
    filtered_df = calculate_zygosity_confidence(edits_df)  # Assumes this function updates the DataFrame in-place
    zygosity_confidence = filtered_df['Class_Conf'].mean()  # Average confidence across all entries


    return {
        "Class_Conf": zygosity_confidence,
        "peaks": ','.join([str(peak) for peak in peak_proportions]),
        "edition_peak_count": len(significant_peaks),
        "max_peak": max_peak,
        "av_peak": np.mean(peak_proportions) if peak_proportions else 0,
        "peak_occupancy": peak_occupancy,
        "clonality": clonality
    }



def parse_indels(csv_path):
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading the CSV file: {e}")
        sys.exit(1)

    # Ensure string type for columns that will use string methods
    for column in ['pre_ins_nt', 'ins_nt', 'post_ins_nt']:
        df[column] = df[column].astype(str)

    # Processing insertions: filter out 'N' and check if DataFrame is empty
    ins_df = df[df['Modification'] == 'ins']
    ins_df = ins_df[
        ~(ins_df['pre_ins_nt'].str.contains('N') |
          ins_df['ins_nt'].str.contains('N') |
          ins_df['post_ins_nt'].str.contains('N'))
    ]

    if ins_df.empty:
        grouped_ins = pd.DataFrame(columns=['Start', 'Length', 'pre_ins_nt', 'ins_nt', 'post_ins_nt', 'Count'])
    else:
        grouped_ins = ins_df.groupby(['Start', 'Length', 'pre_ins_nt', 'ins_nt', 'post_ins_nt']).size().reset_index(name='Count')

    # Process deletions: Filter by 'del'/'delin' and handle empty DataFrame
    dels_df = df[df['Modification'].isin(['del', 'delin'])]
    if dels_df.empty:
        grouped_dels = pd.DataFrame(columns=['Start', 'Length', 'Count'])
    else:
        grouped_dels = dels_df.groupby(['Start', 'Length']).size().reset_index(name='Count')

    return grouped_dels, grouped_ins


def additional_indels_cols(df):
    # Calculate percentages for in-frame and out-of-frame deletions and insertions
    # Initialize the columns to store the sums of outframe and inframe deletions and insertions
    df['Outframe'] = 0
    df['Inframe'] = 0

    # Check if the necessary base columns exist before attempting calculations
    required_columns = ['Dels_inframe', 'Dels_outframe', 'Ins_inframe', 'Ins_outframe', 'Total Reads']
    if all(col in df.columns for col in required_columns):
        # Aggregate inframe and outframe mutations
        df['Inframe'] = df['Dels_inframe'] + df['Ins_inframe']
        df['Outframe'] = df['Dels_outframe'] + df['Ins_outframe']

        # Calculate the percentage for Inframe and Outframe
        df['% Inframe'] = (df['Inframe'] / df['Total Reads']).fillna(0) * 100
        df['% Outframe'] = (df['Outframe'] / df['Total Reads']).fillna(0) * 100

        # Handle any potential division by zero issues by replacing infinities with zero
        df['% Inframe'] = df['% Inframe'].replace([np.inf, -np.inf], 0)
        df['% Outframe'] = df['% Outframe'].replace([np.inf, -np.inf], 0)
    else:
        # If any essential columns are missing, set default percentage values to zero
        df['% Inframe'] = 0
        df['% Outframe'] = 0

    # Now, df contains two new columns: '% Inframe' and '% Outframe' with the calculated percentages.
    return df



def main():

    min_read_threshold = 200


    indel_csv_path = $indels
    edits_csv_path = $edition

    grouped_dels, grouped_ins = parse_indels(indel_csv_path)  
    # Load edits data
    edits_df = pd.read_csv(edits_csv_path)
    # Rename the first column which currently has a blank name
    edits_df.rename(columns={edits_df.columns[0]: 'Sample'}, inplace=True)
    edits_df = parse_edits_csv(edits_df)
    edits_df = additional_indels_cols(edits_df)
    # Initialise zero values in new columns
    edits_df = edits_df.assign(
        Class_Conf=0,
        max_peak=0,
        av_peak=0,
        peak_occupancy=0
    )

    analysis_results = analyze_clonality(grouped_dels, grouped_ins, edits_df, min_read_threshold)
    # Combine with analysis results
    for key in analysis_results:
        edits_df[key] = analysis_results[key]

    outfile = edits_csv_path.replace('.csv','_classified.csv')
    edits_df.to_csv(outfile)
    print(edits_df)


if __name__ == "__main__":
    main()
