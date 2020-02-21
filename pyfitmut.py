import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import minimize
from scipy.optimize import Bounds
import math
from tqdm import tqdm
import argparse


read_num_t1_global = None
read_depth_seq_t1t2_global = None
cell_depth_seq_t1t2_global = None
t_seq_t1t2_global = None
PdfConditional_measure_global = None
read_num_t2_right_global = None

read_num_measure_global = None
cell_num_measure_global = None
read_depth_seq_global = None
cell_depth_seq_global = None
x_mean_global = None
t_seq_global = None
seq_num_global = None
kappa_global = None
c_global = None
Ub_global = None


def fun_filter(read_num_seq):
    # ------------------------------------------------------------------------------------------------------------------
    # A SUB-FUNCTION CALLED BY MAIN FUNCTION main() TO REMOVE NOISY LINEAGES THAT WITH SMALL READ NUMBER)
    #
    # INPUTS
    # --read_num_seq: read number per genotype at each sequencing time-point
    #
    # OUTPUTS
    # --fun_filter_output: label (row number) of filtered lineages,
    #                      read number of filtered lineages
    # ------------------------------------------------------------------------------------------------------------------
    seq_num = np.shape(read_num_seq)[1]
    read_num_seq_DataFrame = pd.DataFrame(read_num_seq)
    read_num_seq_DataFrame = read_num_seq_DataFrame.sort_values(by=list(np.arange(seq_num - 1, -1, -1)),
                                                                ascending=False)
    lineage_id = np.array(read_num_seq_DataFrame.index)
    read_num_seq_DataFrame = read_num_seq_DataFrame.reset_index(drop=True)
    read_num_seq = np.array(read_num_seq_DataFrame, dtype=float)

    # Delete lineages with small read number
    # pos_part = [i for i in range(np.shape(read_num_seq)[0]) for k in range(np.shape(read_num_seq)[1]-1)
    #             if read_num_seq[i, k+1]/(read_num_seq[i, k]+1) >= 50]
    pos = np.logical_not((np.sum(read_num_seq >= 5, axis=1) < 5) + (np.max(read_num_seq, axis=1) < 10)
                         + (np.sum(read_num_seq == 0, axis=1) >= 5) * (read_num_seq[:, -1] >= 1e4))
    read_num_seq = read_num_seq[pos]
    lineage_id = lineage_id[pos]

    filter_output = {'Lineage_ID': lineage_id, 'Read_Number': read_num_seq}

    return filter_output


def fun_mean_fitness(x):
    # ------------------------------------------------------------------------------------------------------------------
    # A SUB-FUNCTION CALLED BY MAIN FUNCTION main() TO ESTIMATE THE MEAN FITNESS AND KAPPA VALUE OF THE POPULATION
    #
    # INPUTS
    # --x: the kappa value and and the mean fitness, [kappa, x_mean]
    #
    # OUTPUTS
    # --square_sum_pdf: the square sum of the difference between theoretical and measured conditional distributions
    # ------------------------------------------------------------------------------------------------------------------
    global read_num_t1_global
    global read_depth_seq_t1t2_global
    global cell_depth_seq_t1t2_global
    global t_seq_t1t2_global
    global PdfConditional_measure_global
    global read_num_t2_right_global

    kappa = x[0]
    x_mean = x[1]

    read_num_t2_est = (read_num_t1_global * np.exp(-x_mean * (t_seq_t1t2_global[1] - t_seq_t1t2_global[0])) *
                       cell_depth_seq_t1t2_global[0] * read_depth_seq_t1t2_global[1] /
                       (read_depth_seq_t1t2_global[0] * cell_depth_seq_t1t2_global[1]))

    read_num_t2 = np.arange(1, read_num_t2_right_global + 0.001)

    PdfConditional_theory = np.multiply(np.sqrt(np.sqrt(read_num_t2_est)
                                                / (4 * np.pi * kappa * np.power(read_num_t2, 1.5))),
                                        np.exp(-np.power(np.sqrt(read_num_t2_est)
                                                         - np.sqrt(read_num_t2), 2) / kappa))

    square_sum_pdf = np.sum(np.power(PdfConditional_theory - PdfConditional_measure_global, 2))*1e-7

    return square_sum_pdf


def fun_likelihood_lineage_mut(x):
    # ------------------------------------------------------------------------------------------------------------------
    # A SUB-FUNCTION CALLED BY MAIN FUNCTION main() TO CALCULATE THE SUM OF THE NEGATIVE LOG LIKELIHOOD VALUE OF EACH
    # LINEAGE (might need more explanation ...)
    #
    # INPUTS
    # --x: the fitness and the establishment time of the mutation, [s, tau]
    #
    # OUTPUTS
    # --likelihood_log: the negative log likelihood value of the lineage
    # ------------------------------------------------------------------------------------------------------------------

    global read_num_measure_global  # measured number of reads
    global cell_num_measure_global  # measured number of cells
    global read_depth_seq_global
    global cell_depth_seq_global
    global x_mean_global
    global t_seq_global
    global seq_num_global
    global kappa_global
    global c_global
    global Ub_global

    s = x[0]
    tau = x[1]

    # --cell_num_theory_neutral: estimated cell number (non-mutant cells) of a neutral lineage
    cell_num_theory_neutral = np.zeros(seq_num_global, dtype=float)
    cell_num_theory_neutral[0] = cell_num_measure_global[0]
    for k in range(1, seq_num_global):
        cell_num_theory_neutral[k] = np.multiply(cell_num_theory_neutral[k - 1],
                                                 np.exp(-(x_mean_global[k] + x_mean_global[k - 1])
                                                        * (t_seq_global[k] - t_seq_global[k - 1]) / 2))

    # --cell_num_theory_adaptive: estimated cell number (both non-mutant ans mutant cells) of a adaptive lineage
    if tau > 0:
        pos_index = [[k, k + 1] for k, ele in enumerate(t_seq_global)
                     if t_seq_global[k] < tau <= t_seq_global[k + 1]][0]
        x_mean_tau = np.interp(tau, [t_seq_global[pos_index[0]], t_seq_global[pos_index[1]]],
                               [x_mean_global[pos_index[0]], x_mean_global[pos_index[1]]])
        established_size_cell_num = (np.true_divide(c_global, s - x_mean_tau) * ((s - x_mean_tau) >= 0.005)
                                     + (c_global / 0.005) * ((s - x_mean_tau) < 0.005))

        cell_num_theory_adaptive_mutant = np.zeros(seq_num_global, dtype=float)
        cell_num_theory_adaptive_mutant[pos_index[1]] = np.multiply(established_size_cell_num,
                                                                    np.exp(s * (t_seq_global[pos_index[1]] - tau)
                                                                           - (x_mean_global[pos_index[1]] + x_mean_tau)
                                                                           * (t_seq_global[pos_index[1]] - tau) / 2))
        if pos_index[1] + 1 < seq_num_global:
            for k in range(pos_index[1] + 1, seq_num_global):
                cell_num_theory_adaptive_mutant[k] = np.multiply(cell_num_theory_adaptive_mutant[k - 1],
                                                                 np.exp(s * (t_seq_global[k] - t_seq_global[k - 1])
                                                                        - (x_mean_global[k] + x_mean_global[k - 1])
                                                                        * (t_seq_global[k] - t_seq_global[k - 1]) / 2))
    elif tau <= 0:
        x_mean_tau = x_mean_global[0]
        established_size_cell_num = (np.true_divide(c_global, s - x_mean_tau) * ((s - x_mean_tau) >= 0.005)
                                     + (c_global / 0.005) * ((s - x_mean_tau) < 0.005))
        cell_num_theory_adaptive_mutant = np.zeros(seq_num_global, dtype=float)
        cell_num_theory_adaptive_mutant[0] = np.multiply(established_size_cell_num,
                                                         np.exp(s * (t_seq_global[0] - tau)
                                                                - (x_mean_global[0] + x_mean_tau)
                                                                * (t_seq_global[0] - tau) / 2))
        for k in range(1, seq_num_global):
            cell_num_theory_adaptive_mutant[k] = np.multiply(cell_num_theory_adaptive_mutant[k - 1],
                                                             np.exp(s * (t_seq_global[k] - t_seq_global[k - 1])
                                                                    - (x_mean_global[k] + x_mean_global[k - 1])
                                                                    * (t_seq_global[k] - t_seq_global[k - 1]) / 2))

    cell_num_theory_adaptive = np.zeros(seq_num_global, dtype=float)
    cell_num_theory_adaptive[0] = cell_num_measure_global[0]
    for k in range(1, seq_num_global):
        tempt_mutant = np.min([cell_num_theory_adaptive[k - 1], cell_num_theory_adaptive_mutant[k - 1]])
        tempt_nonmutant = cell_num_theory_adaptive[k - 1] - tempt_mutant
        cell_num_theory_adaptive[k] = (np.multiply(tempt_mutant, np.exp(s * (t_seq_global[k] - t_seq_global[k - 1])
                                                                        - (x_mean_global[k] + x_mean_global[k - 1])
                                                                        * (t_seq_global[k] - t_seq_global[k - 1]) / 2))
                                       + np.multiply(tempt_nonmutant,
                                                     np.exp(-(x_mean_global[k] + x_mean_global[k - 1])
                                                            * (t_seq_global[k] - t_seq_global[k - 1]) / 2)))

    ratio = np.true_divide(read_depth_seq_global, cell_depth_seq_global)
    read_num_theory_neutral = np.multiply(cell_num_theory_neutral, ratio)
    read_num_theory_adaptive = np.multiply(cell_num_theory_adaptive, ratio)

    # Calculate likelihood (in log)
    delta_s = 0.005
    tempt1 = Ub_global * np.exp(-s / 0.1) / 0.1 * delta_s
    probability_prior_adaptive_log = (np.log(tempt1) + np.log(s / sp.special.gamma(read_num_measure_global[0] * tempt1))
                                      - read_num_measure_global[0] * tempt1 * s / c_global * tau - np.exp(
                -s * tau))  # need change
    probability_prior_neutral_log = 0
    likelihood_adaptive_log = (1 / 4 * np.log(read_num_theory_adaptive)
                               - np.power(np.sqrt(read_num_measure_global)
                                          - np.sqrt(read_num_theory_adaptive), 2) / kappa_global)
    likelihood_neutral_log = (1 / 4 * np.log(read_num_theory_neutral)
                              - np.power(np.sqrt(read_num_measure_global)
                                         - np.sqrt(read_num_theory_neutral), 2) / kappa_global)

    likelihood_log = (np.sum(likelihood_adaptive_log - likelihood_neutral_log)
                      + probability_prior_adaptive_log - probability_prior_neutral_log)

    return -likelihood_log


def main():
    # ------------------------------------------------------------------------------------------------------------------
    # ESTIMATE FITNESS AND ESTABLISHMENT TIME OF EACH SPONTANEOUS ADAPTIVE MUTATION IN COMPETITIVE POOLED GROWTH OF A
    # ISOGENIC POPULATION
    #
    # OPTIONS
    # --input: a .csv file, with each column being the read number per genotype at each sequenced time-point
    # --t_seq: sequenced time-points in number of generations (format: 0 t1 t2 ...)
    # --output_filename: prefix of output .csv files (default: output)
    #
    # OUTPUTS
    # output_filename_FitMut_Result.csv: 1st column: estimated fitness of each genotype, [x1, x2, ...],
    #                                    2nd column: log likelihood value of each genotype, [f1, f2, ...],
    #                                    3rd column: estimated mean fitness per sequenced time-point
    #                                                [x_mean(0), x_mean(t1), ...],
    #                                    4th column+: estimated reads number per genotype per sequencing time-point,
    #                                                 with each time-point being a column
    # ------------------------------------------------------------------------------------------------------------------
    global read_num_t1_global
    global read_depth_seq_t1t2_global
    global cell_depth_seq_t1t2_global
    global t_seq_t1t2_global
    global PdfConditional_measure_global
    global read_num_t2_right_global

    global read_num_measure_global  # measured number of reads
    global cell_num_measure_global  # measured number of cells
    global read_depth_seq_global
    global cell_depth_seq_global
    global x_mean_global
    global t_seq_global
    global seq_num_global
    global kappa_global
    global c_global
    global Ub_global

    parser = argparse.ArgumentParser(description='Estimate fitness and establishment time of each spontanuous adaptive '
                                                 'mutations in a competitive pooled growth experiment',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', type=str, help='a .csv file: with each column being the read number per '
                                                        'genotype at each sequenced time-point')
    parser.add_argument('-t', '--t_seq', nargs='*', type=float, help='sequenced time-points in number of generations')
    parser.add_argument('-n', '--cell_depth_seq', nargs='*', type=float,
                        help='number of cells of the population for each time-point')
    parser.add_argument('-u', '--mutation_rate', type=float, default=1e-5, help='total beneficial mutation rate')
    parser.add_argument('-c', '--c', type=float, default=0.5, help='a noise parameter that characterizes the cell '
                                                                   'growth and cell transfer')  # might need change

    parser.add_argument('-o', '--output_filename', type=str, default='output', help='prefix of output .csv files')

    args = parser.parse_args()
    read_num_seq = np.array(pd.read_csv(args.input, header=None), dtype=float)
    t_seq = np.array(args.t_seq, dtype=float)
    cell_depth_seq = np.array(args.cell_depth_seq, dtype=float)
    Ub = args.mutation_rate
    c = args.c
    output_filename = args.output_filename

    lineages_num_unfiltered = np.shape(read_num_seq)[0]

    # remove very noisy lineages
    filter_output = fun_filter(read_num_seq)

    # estimate mean fitness and kappa value
    read_num_t1_left = 20
    read_num_t1_right = 30
    opt_ini = [2.5, 0]

    read_num_seq = filter_output['Read_Number']
    read_num_seq[read_num_seq == 0] = 1e-7
    read_depth_seq = np.sum(read_num_seq, axis=0)
    cell_depth_seq = np.array(cell_depth_seq, dtype=float)

    read_num_t2_left = 0
    read_num_t2_right = 4 * read_num_t1_right
    x0 = opt_ini
    seq_num = len(t_seq)
    kappa_seq = np.zeros((read_num_t1_right - read_num_t1_left, seq_num - 1), dtype=float)
    x_mean_seq = np.zeros((read_num_t1_right - read_num_t1_left, seq_num - 1), dtype=float)

    for read_num_t1 in range(read_num_t1_left, read_num_t1_right):
        read_num_t1_global = read_num_t1
        for k in range(seq_num - 1):
            read_depth_seq_t1t2_global = read_depth_seq[k:k + 2]
            cell_depth_seq_t1t2_global = cell_depth_seq[k:k + 2]
            t_seq_t1t2_global = t_seq[k:k + 2]
            read_num_t2_right_global = read_num_t2_right
            pos = read_num_seq[:, k] == read_num_t1
            PdfConditional_measure_global = np.histogram(read_num_seq[pos, k + 1],
                                                         bins=np.arange(read_num_t2_left, read_num_t2_right + 0.001),
                                                         density=True)[0]

            opt_output = minimize(fun_mean_fitness, x0, method='BFGS',
                                  options={'disp': False, 'maxiter': 1000, 'gtol': 1e-10, 'eps': 1e-10})
            kappa_seq[read_num_t1 - read_num_t1_left, k], x_mean_seq[read_num_t1 - read_num_t1_left, k] = opt_output.x

    mean_fitness_output = {'Mean_Fitness': [np.max((i, 0)) for i in np.mean(x_mean_seq, axis=0)],
                           'Kappa': np.mean(kappa_seq, axis=0)}

    # estimate fitness and establishment time of each adaptive mutation
    t_seq = t_seq[:-1]
    # seq_num = len(t_seq)
    cell_depth_seq = cell_depth_seq[:-1]
    read_num_seq = read_num_seq[:, :-1]
    x_mean = mean_fitness_output['Mean_Fitness']
    kappa = mean_fitness_output['Kappa']

    read_num_seq = np.array(read_num_seq, dtype=float)
    read_num_seq[read_num_seq == 0] = 1e-7
    cell_depth_seq = np.array(cell_depth_seq, dtype=float)

    read_depth_seq_global = np.sum(read_num_seq, axis=0)
    cell_depth_seq_global = cell_depth_seq
    x_mean_global = x_mean
    t_seq_global = t_seq
    seq_num_global = len(t_seq)
    kappa_global = kappa
    c_global = c
    Ub_global = Ub

    cell_num_seq = read_num_seq / read_depth_seq_global * cell_depth_seq
    lineages_num = np.shape(read_num_seq)[0]
    result_output = {'Mutation_Fitness': np.zeros(lineages_num_unfiltered, dtype=float),
                     'Establishment_Time': np.zeros(lineages_num_unfiltered, dtype=float),
                     'Likelihood_Log': np.zeros(lineages_num_unfiltered, dtype=float),
                     'Mean_Fitness': x_mean}

    lineage_id = filter_output['Lineage_ID']

    x0 = [0.005, 1]
    bounds = Bounds([0.001, -150], [0.5, math.floor(t_seq[-1] - 1)])
    for i in tqdm(range(int(lineages_num / 5))):
        # if i != 0 and i % 1000 == 0:
        #     print(str(i) + ' lineages parsed...')
        read_num_measure_global = read_num_seq[i, :]
        cell_num_measure_global = cell_num_seq[i, :]
        opt_output = minimize(fun_likelihood_lineage_mut, x0, method='L-BFGS-B', bounds=bounds, tol=1e-10)
        result_output['Mutation_Fitness'][lineage_id[i]] = opt_output.x[0]
        result_output['Establishment_Time'][lineage_id[i]] = opt_output.x[1]
        result_output['Likelihood_Log'][lineage_id[i]] = -opt_output.fun

    result_output_csv = pd.DataFrame(result_output)
    result_output_csv.to_csv(output_filename + '_MutSeq_Result.csv', index=False)


if __name__ == "__main__":
    main()
