import numpy as np
import pandas as pd
import copy
import itertools
import csv
import argparse
from tqdm import tqdm


def random_arbitrary(n, bins_edge, freq_bin):
    # generate n random numbers that follow an arbitrary probability distribution
    # the distribution is given by bins_edge = [x0, x1, x2, ..., xn]
    # and counts frequency in each bin freq_bin = [f0, f1, ..., f_{n-1}] f0 in [x0, x1)
    freq_bin_accum = np.cumsum(freq_bin)
    n = int(n)
    output = np.zeros(n)
    k = 0
    while k < n:
        tempt = np.random.rand(1)
        pos = np.where(tempt < freq_bin_accum)[0][0]
        l1 = bins_edge[pos]
        l2 = bins_edge[pos + 1]
        output[k] = l1 + (l2 - l1) * tempt
        k += 1
    return output


def main():
    # ------------------------------------------------------------------------------------------------------------------
    # SIMULATED COMPETITIVE POOLED GROWTH OF A ISOGENIC POPULATION WITH SPONTANEOUS MUTATIONS.
    # THESE SIMULATIONS INCLUDE EXPERIMENTAL NOISE SUCH AS GROWTH NOISE, SAMPLING DURING BOTTLENECKS, DNA EXTRACTION,
    # PCR, AND SAMPLING ON SEQUENCER.
    #
    # OPTIONS
    # --input: a .csv file, with the 1st column being initial cell number of each genotype at generation 0, [n1,n2,...],
    #           and the 2nd and 3rd columns defining an probability distribution for the fitness of mutations (the 2nd
    #           column is the bins edges, [x0, x1, x2, ..., xn], and the 3rd columns is the counts frequency fi in each
    #           bin [xi, x_{i+1}) , [f0, f1, ..., f_{n-1}])
    # --t_seq: time-points evaluated in number of generations (format: 0 t1 t2 ...)
    # --read_num_average_seq: average number of reads per genotype for each time-point (format: 0 r1 r2 ...)
    # --noise_option: which types of noise to include in the simulation, default is all sources of noise
    #                 (`default: growth bottleneck_transfer DNA_extraction PCR sequencing`)
    # --dna_copies: average genome copy number per genotype used as template in PCR (default: 500)
    # --pcr_cycles: number of cycles of PCR (default: 25)
    # --maximum_mutation_number: number of maximum mutations allowed for each lineage
    # --mutation_rate: total beneficial mutation rate
    # --output_filename: prefix of output .csv files (default: output)
    #
    # OUTPUTS
    # output_filename_EvoSimulation_Read_Number.csv: read number per genotype for each time-point
    # output_filename_EvoSimulation_Other_Info.csv: a record of all inputs and other simulated outputs
    # ------------------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Simulated competitive pooled growth of a population of genotypes '
                                                 'with different fitnesses',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str,
                        help='a .csv file: 1st column is initial cell number of each genotype at generation 0, 2nd '
                             'column is bin edges of the arbitrary probability distribution for fitness effect of '
                             'mutations, 3rd column is counts frequency in each bin of the arbitrary probability '
                             'distribution for fitness effect of mutations')
    parser.add_argument('-t', '--t_seq', nargs='*', type=int, help='time-points evaluated in number of generations')
    parser.add_argument('-r', '--read_num_average_seq', nargs='*', type=int,
                        help='average number of reads per genotype for each time-point')
    parser.add_argument('-n', '--noise_option', nargs='*', type=str,
                        default=['growth', 'bottleneck_transfer', 'DNA_extraction', 'PCR', 'sequencing'],
                        help='which types of noise to include in the simulation, default is all sources of noise')
    parser.add_argument('-d', '--dna_copies', type=int, default=500,
                        help='average genome copy number per genotype used as template in PCR')
    parser.add_argument('-p', '--pcr_cycles', type=int, default=25, help='number of cycles of PCR')
    parser.add_argument('-m', '--maximum_mutation_number', type=int, default=1,
                        help='number of maximum mutations allowed for each lineage')
    parser.add_argument('-u', '--mutation_rate', type=float, default=1e-5, help='total beneficial mutation rate')
    parser.add_argument('-o', '--output_filename', type=str, default='output', help='prefix of output .csv files')

    args = parser.parse_args()
    csv_input = pd.read_csv(args.input, header=None, keep_default_na=False)
    cell_num_ini = np.array([i for i in list(csv_input[0]) if i != ''], dtype=int)
    bins_edge = np.array([i for i in list(csv_input[1]) if i != ''], dtype=float)
    freq_bin = np.array([i for i in list(csv_input[2]) if i != ''], dtype=float)
    t_seq = np.array(args.t_seq, dtype=int)
    read_num_average_seq = np.array(args.read_num_average_seq, dtype=int)
    noise_option = args.noise_option
    dna_copies = args.dna_copies
    pcr_cycles = args.pcr_cycles
    lineage_maximum_mutation_number = args.maximum_mutation_number
    Ub = args.mutation_rate
    output_filename = args.output_filename

    delta_t = t_seq[1] - t_seq[0]
    lineages_num = np.max(cell_num_ini.shape)
    seq_num = np.max(t_seq.shape)
    evo_num = t_seq[-1]
    c = 0.5

    if ('growth' in noise_option) and ('bottleneck_transfer' in noise_option):
        f1 = np.random.poisson
        f2 = np.random.poisson
    elif ('growth' in noise_option) and ('bottleneck_transfer' not in noise_option):
        f1 = np.random.poisson
        f2 = np.round
    elif ('growth' not in noise_option) and ('bottleneck_transfer' in noise_option):
        f1 = np.round
        f2 = np.random.poisson
    else:
        f1 = np.round
        f2 = np.round

    if ('DNA_extraction' in noise_option) and ('PCR' in noise_option) and ('sequencing' in noise_option):
        f3 = np.random.poisson
        f4 = np.random.poisson
        f5 = np.random.poisson
    elif ('DNA_extraction' in noise_option) and ('PCR' in noise_option) and ('sequencing' not in noise_option):
        f3 = np.random.poisson
        f4 = np.random.poisson
        f5 = np.round
    elif ('DNA_extraction' in noise_option) and ('PCR' not in noise_option) and ('sequencing' in noise_option):
        f3 = np.random.poisson
        f4 = np.round
        f5 = np.random.poisson
    elif ('DNA_extraction' not in noise_option) and ('PCR' in noise_option) and ('sequencing' in noise_option):
        f3 = np.round
        f4 = np.random.poisson
        f5 = np.random.poisson
    elif ('DNA_extraction' in noise_option) and ('PCR' not in noise_option) and ('sequencing' not in noise_option):
        f3 = np.random.poisson
        f4 = np.round
        f5 = np.round
    elif ('DNA_extraction' not in noise_option) and ('PCR' in noise_option) and ('sequencing' not in noise_option):
        f3 = np.round
        f4 = np.random.poisson
        f5 = np.round
    elif ('DNA_extraction' not in noise_option) and ('PCR' not in noise_option) and ('sequencing' in noise_option):
        f3 = np.round
        f4 = np.round
        f5 = np.random.poisson
    else:
        f3 = np.round
        f4 = np.round
        f5 = np.round

    x_mean = np.zeros(evo_num + 1)

    if lineage_maximum_mutation_number == 1:  # [s, n, occurrence of mutation, establishment of mutation]
        # Growth phase, with two possible noise involved: cell growth noise, bottleneck cell transfer noise
        cell_num_tempt = {i: {'mut0': cell_num_ini[i], 'mut1': []} for i in range(lineages_num)}
        cell_num_bottleneck_seq = dict()
        cell_num_bottleneck_seq[0] = {i: {'mut0': cell_num_ini[i] * 2 ** delta_t, 'mut1': []} for i in
                                      range(lineages_num)}
        for j in range(1, seq_num):
            cell_num_bottleneck_seq[j] = {i: {'mut0': 0, 'mut1': []} for i in range(lineages_num)}

        for j in tqdm(range(1, evo_num + 1)):
            x_rela0 = np.max([1 / (1 + x_mean[j - 1]), 0])
            for i in range(lineages_num):
                tempt = copy.deepcopy(cell_num_tempt[i])
                tempt1 = tempt['mut1']
                if tempt1:
                    for k in range(len(tempt1)):
                        if tempt1[k][1] > 0:
                            x_rela = np.max([(1 + tempt1[k][0]) / (1 + x_mean[j - 1]), 0])
                            cell_num_tempt[i]['mut1'][k][1] = f1(2 * x_rela * tempt1[k][1])

                tempt0 = cell_num_tempt[i]['mut0']
                if tempt0 > 0:
                    cell_num_tempt[i]['mut0'] = f1(2 * x_rela0 * tempt0)
                    mut1_num = np.random.binomial(tempt0, Ub)
                    if mut1_num > 0:
                        for m in range(mut1_num):
                            cell_num_tempt[i]['mut1'].append(
                                [random_arbitrary(1, bins_edge, freq_bin)[0], 1, j, -1000])

            # -- cell transfer at the bottleneck
            if np.mod(j, delta_t) == 0:
                ind = int(j / delta_t)
                cell_num_bottleneck_seq[ind] = copy.deepcopy(cell_num_tempt)
                for i in range(lineages_num):
                    cell_num_tempt[i]['mut0'] = f2(cell_num_tempt[i]['mut0'] / (2 ** delta_t))
                    tempt1 = cell_num_tempt[i]['mut1']
                    if tempt1:
                        for k in range(len(tempt1)):
                            cell_num_tempt[i]['mut1'][k][1] = f2(tempt1[k][1] / (2 ** delta_t))

            # -- relabel establishment of mutations
            for i in range(lineages_num):
                tempt1 = cell_num_tempt[i]['mut1']
                if tempt1:
                    for k in range(len(tempt1)):
                        if (tempt1[k][1] > c / (tempt1[k][0] + 1e-10)) and (tempt1[k][3] == -1000):
                            cell_num_tempt[i]['mut1'][k][3] = j
                        elif (tempt1[k][1] < c / (tempt1[k][0] + 1e-10)) and (tempt1[k][3] > -1000):
                            cell_num_tempt[i]['mut1'][k][3] = -1000

            # -- relabel establishment of mutations
            if np.mod(j, delta_t) == 0:
                ind = int(j / delta_t)
                for i in range(lineages_num):
                    tempt1 = cell_num_tempt[i]['mut1']
                    if tempt1:
                        for k in range(len(tempt1)):
                            cell_num_bottleneck_seq[ind][i]['mut1'][k][3] = tempt1[k][3]

            # -- mean fitness
            cell_num_total_mut1 = 0
            mean_tempt_mut1 = 0
            cell_num_total_mut0 = sum([cell_num_tempt[i]['mut0'] for i in range(lineages_num)])
            for i in range(lineages_num):
                tempt1 = cell_num_tempt[i]['mut1']
                if tempt1:
                    for k in range(len(tempt1)):
                        cell_num_total_mut1 += tempt1[k][1]
                        mean_tempt_mut1 += tempt1[k][0] * tempt1[k][1]
            x_mean[j] = mean_tempt_mut1 / (cell_num_total_mut0 + cell_num_total_mut1)

        # After-growth phase, with three possible noise involved: DNA extraction noise, PCR noise, sequencing noise
        cell_num_seq_array = np.zeros((lineages_num, seq_num))
        for j in range(seq_num):
            for i in range(lineages_num):
                tempt = cell_num_bottleneck_seq[j][i]
                tempt0 = tempt['mut0']
                tempt1 = [tempt['mut1'][k][1] for k in range(len(tempt['mut1']))]
                cell_num_seq_array[i, j] = tempt0 + sum(tempt1)

        dna_num_seq_array = f3(np.true_divide(cell_num_seq_array, np.sum(cell_num_seq_array, axis=0))
                               * dna_copies * lineages_num)
        pcr_num_seq_array = dna_num_seq_array
        for i in range(pcr_cycles):
            pcr_num_seq_array = f4(2 * pcr_num_seq_array)
        read_num_seq_array = np.multiply(np.true_divide(pcr_num_seq_array, np.sum(pcr_num_seq_array, axis=0)),
                                         read_num_average_seq * lineages_num)
        read_num_seq_array = f5(read_num_seq_array)

        tempt = cell_num_bottleneck_seq[seq_num - 1]
        mut1_info = [[tempt[i]['mut1'][k][0], tempt[i]['mut1'][k][2], tempt[i]['mut1'][k][3], i, k]
                     for i in range(lineages_num) for k in range(len(tempt[i]['mut1']))
                     if tempt[i]['mut1'][k][3] > -1000 and tempt[i]['mut1'][k][1] > 0]

        evo_simulator_output = {'Read_Number': read_num_seq_array,
                                'Other_Info': {'Mutation1_Information': mut1_info,
                                               'Mean_Fitness': x_mean[t_seq],
                                               'Time_Points': t_seq,
                                               'Average_Read_Depth': read_num_average_seq,
                                               'Noise': noise_option,
                                               'gDNA_Copies': [dna_copies],
                                               'PCR_cycles': [pcr_cycles],
                                               'Lineage_Maximum_Mut_Number': [lineage_maximum_mutation_number],
                                               'Ub': [Ub]}}

    elif lineage_maximum_mutation_number == 2:   # [label, s, n, occurrence of mutation, establishment of mutation]
        # Growth phase, with two possible noise involved: cell growth noise, bottleneck cell transfer noise
        cell_num_tempt = {i: {'mut0': cell_num_ini[i], 'mut1': [], 'mut2': []} for i in range(lineages_num)}
        cell_num_bottleneck_seq = dict()
        cell_num_bottleneck_seq[0] = {i: {'mut0': cell_num_ini[i] * 2 ** delta_t, 'mut1': [], 'mut2': []} for i in
                                      range(lineages_num)}
        for j in range(1, seq_num):
            cell_num_bottleneck_seq[j] = {i: {'mut0': 0, 'mut1': [], 'mut2': []} for i in range(lineages_num)}

        for j in tqdm(range(1, evo_num + 1)):
            x_rela0 = np.max([1 / (1 + x_mean[j - 1]), 0])
            for i in range(lineages_num):
                tempt = copy.deepcopy(cell_num_tempt[i])
                tempt1 = tempt['mut1']
                tempt2 = tempt['mut2']
                if tempt2:
                    for k in range(len(tempt2)):
                        if tempt2[k][2] > 0:
                            x_rela = np.max([(1 + tempt1[tempt2[k][0]][0] + tempt2[k][1]) / (1 + x_mean[j - 1]), 0])
                            cell_num_tempt[i]['mut2'][k][2] = f1(2 * x_rela * tempt2[k][2])

                if tempt1:
                    for k in range(len(tempt1)):
                        if tempt1[k][1] > 0:
                            x_rela = np.max([(1 + tempt1[k][0]) / (1 + x_mean[j - 1]), 0])
                            cell_num_tempt[i]['mut1'][k][1] = f1(2 * x_rela * tempt1[k][1])
                            mut2_num = np.random.binomial(tempt1[k][1], Ub)
                            if mut2_num > 0:
                                for m in range(mut2_num):
                                    cell_num_tempt[i]['mut2'].append(
                                        [k, random_arbitrary(1, bins_edge, freq_bin)[0], 1, j, -1000])

                tempt0 = cell_num_tempt[i]['mut0']
                if tempt0 > 0:
                    cell_num_tempt[i]['mut0'] = f1(2 * x_rela0 * tempt0)
                    mut1_num = np.random.binomial(tempt0, Ub)
                    if mut1_num > 0:
                        for m in range(mut1_num):
                            cell_num_tempt[i]['mut1'].append(
                                [random_arbitrary(1, bins_edge, freq_bin)[0], 1, j, -1000])

            # -- cell transfer at the bottleneck
            if np.mod(j, delta_t) == 0:
                ind = int(j / delta_t)
                cell_num_bottleneck_seq[ind] = copy.deepcopy(cell_num_tempt)
                for i in range(lineages_num):
                    cell_num_tempt[i]['mut0'] = f2(cell_num_tempt[i]['mut0'] / (2 ** delta_t))
                    tempt1 = cell_num_tempt[i]['mut1']
                    tempt2 = cell_num_tempt[i]['mut2']
                    if tempt1:
                        for k in range(len(tempt1)):
                            cell_num_tempt[i]['mut1'][k][1] = f2(tempt1[k][1] / (2 ** delta_t))
                    if tempt2:
                        for k in range(len(tempt2)):
                            cell_num_tempt[i]['mut2'][k][2] = f2(tempt2[k][2] / (2 ** delta_t))

            # -- relabel establishment of mutations
            for i in range(lineages_num):
                tempt1 = cell_num_tempt[i]['mut1']
                tempt2 = cell_num_tempt[i]['mut2']
                if tempt1:
                    for k in range(len(tempt1)):
                        if (tempt1[k][1] > c / (tempt1[k][0] + 1e-10)) and (tempt1[k][3] == -1000):
                            cell_num_tempt[i]['mut1'][k][3] = j
                        elif (tempt1[k][1] < c / (tempt1[k][0] + 1e-10)) and (tempt1[k][3] > -1000):
                           cell_num_tempt[i]['mut1'][k][3] = -1000
                if tempt2:
                    for k in range(len(tempt2)):
                        if (tempt2[k][2] > c / (tempt1[tempt2[k][0]][0] + tempt2[k][1] + 1e-10)) and (
                                tempt2[k][4] == -1000):
                            cell_num_tempt[i]['mut2'][k][4] = j
                        elif (tempt2[k][2] < c / (tempt1[tempt2[k][0]][0] + tempt2[k][1] + 1e-10)) and (tempt2[k][4] > -1000):
                           cell_num_tempt[i]['mut2'][k][4] = -1000

            # -- relabel establishment of mutations
            if np.mod(j, delta_t) == 0:
                ind = int(j / delta_t)
                for i in range(lineages_num):
                    tempt1 = cell_num_tempt[i]['mut1']
                    tempt2 = cell_num_tempt[i]['mut2']
                    if tempt1:
                        for k in range(len(tempt1)):
                            cell_num_bottleneck_seq[ind][i]['mut1'][k][3] = tempt1[k][3]
                    if tempt2:
                        for k in range(len(tempt2)):
                            cell_num_bottleneck_seq[ind][i]['mut2'][k][4] = tempt2[k][4]

            # -- mean fitness
            cell_num_total_mut1 = 0
            cell_num_total_mut2 = 0
            mean_tempt_mut1 = 0
            mean_tempt_mut2 = 0
            cell_num_total_mut0 = sum([cell_num_tempt[i]['mut0'] for i in range(lineages_num)])
            for i in range(lineages_num):
                tempt1 = cell_num_tempt[i]['mut1']
                tempt2 = cell_num_tempt[i]['mut2']
                if tempt1:
                    for k in range(len(tempt1)):
                        cell_num_total_mut1 += tempt1[k][1]
                        mean_tempt_mut1 += tempt1[k][0] * tempt1[k][1]
                if tempt2:
                    for k in range(len(tempt2)):
                        cell_num_total_mut2 += cell_num_tempt[i]['mut2'][k][2]
                        mean_tempt_mut2 += (tempt1[tempt2[k][0]][0] + tempt2[k][1]) * tempt2[k][2]
            x_mean[j] = (mean_tempt_mut1 + mean_tempt_mut2) / (cell_num_total_mut0 + cell_num_total_mut1
                                                               + cell_num_total_mut2)

        # After-growth phase, with three possible noise involved: DNA extraction noise, PCR noise, sequencing noise
        cell_num_seq_array = np.zeros((lineages_num, seq_num))
        for j in range(seq_num):
            for i in range(lineages_num):
                tempt = cell_num_bottleneck_seq[j][i]
                tempt0 = tempt['mut0']
                tempt1 = [tempt['mut1'][k][1] for k in range(len(tempt['mut1']))]
                tempt2 = [tempt['mut2'][k][2] for k in range(len(tempt['mut2']))]
                cell_num_seq_array[i, j] = tempt0 + sum(tempt1) + sum(tempt2)

        dna_num_seq_array = f3(np.true_divide(cell_num_seq_array, np.sum(cell_num_seq_array, axis=0))
                               * dna_copies * lineages_num)
        pcr_num_seq_array = dna_num_seq_array
        for i in range(pcr_cycles):
            pcr_num_seq_array = f4(2 * pcr_num_seq_array)
        read_num_seq_array = np.multiply(np.true_divide(pcr_num_seq_array, np.sum(pcr_num_seq_array, axis=0)),
                                         read_num_average_seq * lineages_num)
        read_num_seq_array = f5(read_num_seq_array)

        tempt = cell_num_bottleneck_seq[seq_num - 1]
        mut1_info = [[tempt[i]['mut1'][k][0], tempt[i]['mut1'][k][2], tempt[i]['mut1'][k][3], i, k]
                     for i in range(lineages_num) for k in range(len(tempt[i]['mut1']))
                     if tempt[i]['mut1'][k][3] > -1000 and tempt[i]['mut1'][k][1] > 0]
        mut2_info = [[tempt[i]['mut2'][k][1], tempt[i]['mut2'][k][3], tempt[i]['mut2'][k][4], i,
                      tempt[i]['mut2'][k][0], k] for i in range(lineages_num) for k in range(len(tempt[i]['mut2']))
                     if tempt[i]['mut2'][k][4] > -1000 and tempt[i]['mut2'][k][2] > 0]

        evo_simulator_output = {'Read_Number': read_num_seq_array,
                                'Other_Info': {'Mutation1_Information': mut1_info,
                                               'Mutation2_Information': mut2_info,
                                               'Mean_Fitness': x_mean[t_seq],
                                               'Time_Points': t_seq,
                                               'Average_Read_Depth': read_num_average_seq,
                                               'Noise': noise_option,
                                               'gDNA_Copies': [dna_copies],
                                               'PCR_cycles': [pcr_cycles],
                                               'Lineage_Maximum_Mut_Number': [lineage_maximum_mutation_number],
                                               'Ub': [Ub]}}

    tempt = pd.DataFrame(evo_simulator_output['Read_Number'])
    tempt.to_csv(output_filename + '_EvoSimulation_Read_Number.csv', index=False, header=False)

    tempt = list(itertools.zip_longest(*list(evo_simulator_output['Other_Info'].values())))
    with open(output_filename + '_EvoSimulation_Other_Info.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(evo_simulator_output['Other_Info'].keys())
        w.writerows(tempt)


if __name__ == "__main__":
    main()