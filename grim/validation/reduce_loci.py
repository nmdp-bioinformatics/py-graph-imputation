import argparse
import json


def write_best_prob_genotype(name_gl, res, fout, numOfResult=10):
    sorted_by_value = sorted(res.items(), key=lambda kv: kv[1], reverse=True)

    # write the output to file
    minBestResult = min(numOfResult, len(sorted_by_value))
    for k in range(minBestResult):
        fout.write(name_gl + ',' + str(sorted_by_value[k][0]) + ',' +
                   str(sorted_by_value[k][1]) + ',' + str(k) + '\n')

def convert_res_of_6_to_5(file_in, file_out):
    dict_res = {}
    with open(file_in) as six_file:
        for line in six_file:
            id, gl, prob, rank = line.strip().split(',')
            gl = ('^').join(gl.split('^')[:-1])
            if not id in dict_res:
                dict_res[id] = {}
            if gl in dict_res[id]:
                dict_res[id][gl]+= float(prob)
            else:
                dict_res[id][gl] = float(prob)

    six_file.close()

    f_out = open(file_out, 'w')
    for id in dict_res:
        write_best_prob_genotype(id, dict_res[id], f_out)
    f_out.close()




parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config",
                    required=False,
                    default="../conf/minimal-configuration.json",
                    help="Configuration JSON file",
                    type=str)

args = parser.parse_args()
configuration_file = args.config
with open(configuration_file) as f:
    json_conf = json.load(f)
imputation_out_muug_freqs =  'output/' + json_conf.get("imputation_out_umug_freq_filename")

convert_res_of_6_to_5(imputation_out_muug_freqs, imputation_out_muug_freqs)