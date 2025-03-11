import pandas as pd
import numpy as np
import pathlib as Path

from sklearn.model_selection import RepeatedStratifiedKFold
from scipy.stats import f_oneway,ttest_ind
from statsmodels.stats.multitest import fdrcorrection
from sklearn import metrics

import matplotlib.pyplot as plt

def ROC_curve_analysis(labels:np.ndarray, scores:np.ndarray, title:str):
    """
    Plots the ROC curve using matplotlib
    labes: 1 and 0 indicating positive and negative classes
    scores: continiuos variable the more value is favored for positive class
    """
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
    roc_auc = metrics.auc(fpr, tpr)
    display = metrics.RocCurveDisplay(fpr = fpr,
                                      tpr = tpr,
                                      roc_auc = roc_auc,
                                      estimator_name = title)
    display.plot()
    plt.show()


def intersection(lst1:list, lst2:list):
    """
    The intersection between two lists
    """
    return list(set(lst1) & set(lst2))


def setdiff(lst1:list ,lst2:list):
    """
    The different items in list1 which do not exist in list2
    """
    return(list(set(lst1) - set(lst2)))


def cross_validation_for_one_vs_all_t_test(raw_df:pd.DataFrame,
                                           intensity_file:pd.DataFrame,
                                           meta_data:pd.DataFrame,
                                           Entity:str,
                                           protein_peptides_list,
                                           k_folds:int = 5,
                                           num_repeats:int = 1,
                                           p_value_cutoff:float = .01,
                                           fdr_cutoff:float = .01, 
                                           average_difference_twogroups:float = .75) -> pd.DataFrame:

        
        Y = raw_df['Sarcoma Subtype'].to_numpy()
        X = raw_df.to_numpy()                        

        skf = RepeatedStratifiedKFold(n_splits = k_folds, n_repeats = num_repeats)
        all_data = []
        for i, (train_index, test_index) in enumerate(skf.split(X, Y)):

                # training: calculation of the signatures based on one_vs_all_t_test
                train_set = raw_df.iloc[train_index,:]
                fp_signatures = one_vs_all_t_test(train_set,protein_peptides_list, Entity, 'Sarcoma Subtype')
                fp_signatures['delta'] = abs(fp_signatures['means_group1'] - fp_signatures['means_group2'])
                primary_signatures = fp_signatures[(fp_signatures['p_values'] < p_value_cutoff) &
                                                                (fp_signatures['fdr'] < fdr_cutoff) &
                                                                (fp_signatures['delta'] > average_difference_twogroups)]

                primary_signatures_proteins = list(primary_signatures['Gene Names'])

                # test: calculation of the disease_score
                subset = raw_df.iloc[test_index,:]
                samples = subset['Sample name'].tolist()
                test_patients = intersection(intensity_file.columns,samples)
                results = calculate_entity_score_from_the_signature_proteins(intensity_file, meta_data, primary_signatures_proteins, Entity, test_patients)
                all_data.append(results)

        merged_df = pd.concat(all_data)
        return(merged_df)


def subset_intensity_file_fp(reportDirectort:str,
                             metapath:str,
                             selectedProteins:list,
                             selectedPtientsEntity:str) -> pd.DataFrame:
    """
    Subsets the intensity file based on selected patients and selected proteins
    """
    path_to_file = f'{reportDirectort}/preprocessed_fp.csv'
    intensity_file = pd.read_csv(path_to_file)

    intensity_file = intensity_file[intensity_file['Gene names'].isin(selectedProteins)]
    meta_data = pd.read_excel(metapath)
    selected_patients = list(meta_data.loc[meta_data['Sarcoma Subtype'] == selectedPtientsEntity,'Sample name'])
    selected_df = intensity_file[intersection(intensity_file.columns,selected_patients)]
    selected_df.index = selectedProteins
    return selected_df


def prepareDataframeforTest(reportDirectory:str,
                            fp_pp:str,
                            metapath:str,
                            minimum_patients_per_entity = 20):
    """
    INPUT: - reportDirectory: is the result folder of the wp3 clinical proteomics pipeline
           - fp_pp: fp for full protein and pp for phspho proteins
           - metapath: is the path to the meta file
           - minimum_patients_per_entity: the minimum patients in the entity to consider it as entity

    OUTPUTs: - a dataframe before T_test or ANOVA 
            - a list of the protein names 
    
    """
    path_to_file = f'{reportDirectory}/preprocessed_{fp_pp}.csv'
    intensity_file = pd.read_csv(path_to_file)
    if fp_pp == 'fp':
        intensity_file.index = intensity_file['Gene names']
        protein_peptide = list(intensity_file['Gene names'])
        intensity_file = intensity_file.drop(columns='Gene names')
    else:
        intensity_file.index = intensity_file['Modified sequence']
        protein_peptide = list(intensity_file['Modified sequence'])
        intensity_file = intensity_file.drop(columns='Modified sequence')

    # reading the meta data and considereing the Sarcoma Subtype as the Entity
    meta = pd.read_excel(metapath)
    Entity_count = meta.groupby('Sarcoma Subtype')['Sample name'].count()
    Entities = Entity_count[Entity_count >= minimum_patients_per_entity]     # only the entities with more than n patienst will be used
    Entities = list(Entities.index)
    meta_data = meta[meta['Sarcoma Subtype'].isin(Entities)]
    meta_data = meta_data[['Sample name','Sarcoma Subtype']]

    # Filtering and merging the intensity table to that of the meta data
    patients_for_ANOVA = intersection(intensity_file.columns,meta_data['Sample name'])
    final_inensity = intensity_file[patients_for_ANOVA].transpose()
    final_inensity.columns = intensity_file.index
    ANOVA_df = pd.merge(final_inensity,meta_data,left_index = True,right_on = 'Sample name')
    return ANOVA_df, protein_peptide


def anova_test(inputDF:pd.DataFrame,
               protein_peptide:list,
               metaDataColumn:str) -> pd.DataFrame:
    """
    Performs ANOVA for each protein/peptide row wise 
    INPUT: - a dataframe where the protein/peptides are the rows and
                               the patients are the columns
           - protein_peptide is the list of the columns to do the ANOVA 
           - metaDatacolumns is the column in the inputDf which contains the subtypes for grouping

    """
    # new_Df = new_Df.transpose()
    F_tests,p_values,mean_grps = [],[],[]
    for gene in protein_peptide:
        F,p = 1,1 # a precomputed value for the p and F
        column_names = [gene,metaDataColumn]
        grp = inputDF[column_names]
        grp = grp.dropna()
        #
        result = grp.groupby(metaDataColumn)[gene].apply(list)
        average = grp.groupby(metaDataColumn)[gene].mean()
        average = list(average)
        # more than two groups for each gene
        if (len(grp['localisation'].unique())) > 2:
            F, p = f_oneway(*result)
        F_tests.append(F)
        p_values.append(p)
        mean_grps.append(average)
    p_df = pd.DataFrame(list(zip(F_tests,p_values,mean_grps)))
    p_df.columns = ['F_tests','p_values','means_pergroup']
    p_df['Gene Names'] = protein_peptide
    p_df = p_df[p_df['p_values'].notna()]
    fdr_multi_correction = fdrcorrection(p_df.p_values, alpha=0.05, method='indep', is_sorted=False)
    p_df['fdr'] = list(fdr_multi_correction[1])
    return p_df


def t_test(x:tuple, y:tuple):
    """ Performs t_test """
    t_test = ttest_ind(x, y)
    F = list(t_test)[0]
    p = list(t_test)[1]
    return F,p

def get_final_signatures_from_t_test_results_by_criteria(signatures_df:pd.DataFrame,
                                                         p_value_threshold:float = .01,
                                                         fdr_threshold:float = .01,
                                                         delta_threshold:float = .75) -> pd.DataFrame:
    return(signatures_df[(signatures_df['p_values'] < p_value_threshold) & (signatures_df['fdr'] < fdr_threshold) & (signatures_df['delta'] > delta_threshold) ])
    


def one_vs_all_t_test(inputDF:pd.DataFrame,
                      protein_peptide:list,
                      favoriteentity:str,
                      metaDataColumn:str) -> pd.DataFrame:
    """
    Performs a one vs all T_test per each protein/peptide for the patients with the favorite entity vs all other entities
    to find out the signature protein for that specific entity 

    :inputDF: a dataframe where the patients are the columns and
                               the proteins/peptides are the rows
    :protein_peptide: the list of the columns to do the t_test based on 
    :favoriteentity: the main group for the t_test i.e: chordoma
    :metaDatacolumns: the column in the inputDf which contains the subtypes for grouping

    """
    F_tests, p_values, mean_grp1, mean_grp2, num_grp1, num_other_grps = [],[],[],[],[],[]
    for gene in protein_peptide:
        F,p = 1,1 # a precomputed value for the p and F
        column_names = [gene,metaDataColumn]
        grp = inputDF[column_names]
        grp = grp.dropna()
        df = grp
        g1 = df[(df[metaDataColumn] == favoriteentity)]  # first group for the test
        g2 = df[(df[metaDataColumn] != favoriteentity)]  # second group for the test
        x = tuple(g1[gene])
        y = tuple(g2[gene])
        F, p = t_test(x,y)
        F_tests.append(F)
        p_values.append(p)
        mean_x = g1[gene].mean()
        mean_y = g2[gene].mean()
        num_grp1.append(len(g1[gene]))
        num_other_grps.append(len(g2[gene]))
        mean_grp1.append(mean_x)
        mean_grp2.append(mean_y)
    p_df = pd.DataFrame(list(zip(F_tests,p_values,mean_grp1,mean_grp2,num_grp1,num_other_grps)))
    p_df.columns = ['t_statistics','p_values','means_group1','means_group2','num_samples_groups_interest','num_sample_other_groups']
    p_df['Gene Names'] = protein_peptide
    p_df = p_df[p_df['p_values'].notna()]
    fdr_multi_correction = fdrcorrection(p_df.p_values, alpha=0.05, method='indep', is_sorted=False)
    p_df['fdr'] = list(fdr_multi_correction[1])
    p_df['delta'] = abs(p_df['means_group1'] - p_df['means_group2'])
    return p_df


def calculate_entity_score_from_the_signature_proteins(intensity_file:pd.DataFrame,
                                                        meta_data:pd.DataFrame,
                                                        list_signature_proteins:pd.Series,
                                                        cohort_of_interest:str,
                                                        test_patients:pd.Series) -> pd.DataFrame:
    """
    Takes the result of signatures through one_vs_all t_test and returns a vector of predicted values for the patients IDs.
    The prediction is based on the model:     scores = [intensity] * W
    W is the vector of weights from the t_test and is calculated based on the average of intensities for the patients of interested 
    cohort with respect to other patients it can  be either +1 or -1.
    The intensities NAs will be imputed based on the average of patients for that protein

    :intensity_file: A dataframe with a columns named Gene Names and the other columns as patients, the rows are proteins
    :meta_data: A dataframe with at least "Sarcoma Subtype" column that contains entities including cohort_of_interest
    :signatures_df: list of proteins 
    :cohort_of_interest: is the cohort used to make the model based on that
    :test_patients: the list of patients to be used in prediction

    OUTPUT:
        dataframe of
        - labels a vector of 0 and 1, real values if the patient is with the cohort it will be 1 other wise it will be 0
        - scores the predicted result of [intenstity] * weights

    """
    test_patients = intersection(intensity_file.columns,test_patients)
    intensity_file = intensity_file[intensity_file['Gene names'].isin(list_signature_proteins) ]
    gene_names = intensity_file['Gene names']
    intensity_file = intensity_file[test_patients]
    intensity_file.index = gene_names
    # calculation of the weights
    interested_patients = list(meta_data.loc[meta_data['Sarcoma Subtype'] == cohort_of_interest, 'Sample name'])
    interested_patients = intersection(interested_patients,test_patients)
    weights = [0 for i in range(len(intensity_file.index))]
    print(len(weights))
    for j in range(len(weights)):
        protein = intensity_file.index[j]
        grp = intensity_file[intensity_file.index == protein]
        grp_mean1 = grp[interested_patients].mean(axis =1,skipna=True)
        grp_mean2 = grp.drop(interested_patients , axis=1).mean(axis=1,skipna=True)
        if (grp_mean1[0] > grp_mean2[0]):
            weights[j] = 1    # for up regulation
        else: 
            weights[j] = -1   # for down regulation 

    labels = [1 if i in interested_patients else 0 for i in intensity_file.columns]
    intensity_file = intensity_file.fillna(intensity_file.mean())  # imputation for each gene as the average of the intensities across all patients
    weight_array = pd.DataFrame(weights).to_numpy()
    scores = np.dot(intensity_file.transpose(), weight_array)
    result_df = pd.DataFrame(list(zip(labels,scores)))
    result_df.columns = ['labels','scores']
    result_df.index = intensity_file.columns
    return result_df

if __name__ == '__main__':
    import sys
    import json
    import argparse
    from . import config


    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True,
                        help="Absolute path to configuration file.")

    parser.add_argument("-E", "--entity", required=True,
                    help="the entity to calculate the scores like Chordoma")

    args = parser.parse_args(sys.argv[1:])
    configs = config.load(args.config)
    entity = args.entity
    report_directory = configs['results_folder']
    path_to_file = f'{report_directory}/preprocessed_fp.csv'

    # signature by t_test
    intensity_file = pd.read_csv(path_to_file)
    meta_file_path = configs['metadata_annotation']
    meta_data = pd.read_excel(meta_file_path)
    raw_df , protein_peptides_list = prepareDataframeforTest(report_directory,'fp', meta_file_path, minimum_patients_per_entity = 20)
    fp_signatures = one_vs_all_t_test(raw_df,protein_peptides_list, entity, 'Sarcoma Subtype')
    primary_signatures = get_final_signatures_from_t_test_results_by_criteria(fp_signatures,.01,.01,.75)
    primary_signatures.to_csv(f'{report_directory}/primary_signature_{entity}.csv')
    # Entity scores
    patients_list = intensity_file.filter(regex=r'^[A-Z,a-z]+\d{1,3}-\S+-\S+').columns.tolist()
    entity_scores = calculate_entity_score_from_the_signature_proteins(intensity_file,meta_data,primary_signatures['Gene Names'],entity,patients_list)
    entity_scores.to_csv(f'{report_directory}/entity_scores_{entity}.csv')
    # cross_validation
    merged_df = cross_validation_for_one_vs_all_t_test(raw_df,intensity_file,meta_data,entity,protein_peptides_list)
    ROC_curve_analysis(merged_df['labels'] ,merged_df['scores'],f'ROC_curve_CV_{entity}')

