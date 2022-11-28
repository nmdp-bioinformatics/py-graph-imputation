##Configuration file parameters:

| Parameter | Description |
| --- | --- |
| populations | The population to consider them frequencies. |
| priority | The coefficient values that define the priority matrix. |
| loci_map| Loci full name Mapping for indexes. |
| freq_trim_threshold | The numerator in the frequency threshold. |
| factor_missing_data | factor to haplotype frequency in plan B in missing data case |
| Plan_B_Matrix | matrix arranged by the most probable possibilities for recombination. The first element in the matrix should be the full haplotype. the indexes are corresponding to loci_map|
| planb| True - use plan B anc C. False - use only Plan A. |
| imputation_in_file | Population genotypes input path. |
| imputation_out_umug_freq_filename | imputation output. |
| freq_file | Output file path. |
| number_of_results | A maximum number of corresponding umug/haplotyes pair for each individual. |
| epsilon | Genotypes with smaller probability than this value will be discarded. |
| max_haplotypes_number_in_phase | Limits the number of processed haplotypes in each phase. |
| Plan_A_Matrix | A list of nodes to plan A. The full-locus nodes will be created anyway, whether they are on the list or not. If the list is empty or this field is missing, all nodes and edges will be created. |
| save_space_mode | Reduce options in plan B and C if there are too much alleles. Suitable for 9-locus. Default - False |
