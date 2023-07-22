# Filter Hereditary cancer mutation
# author: Yifei.wan

import csv
import argparse


def check_vista_curation_status(line_object):
    # Check VISta curation status
    vista_curation_status = line_object['VISta_curation_status']

    if vista_curation_status == 'APPROVED':
        return True
    return False


def check_vista_curation_needed(line_object):
    # Check if VISta curation is needed
    vista_curation_status = line_object['VISta_curation_status']

    if vista_curation_status == 'CURATION NEEDED':
        return True
    return False


def check_vista_conclusion_benign(line_object):
    # Check if VISta conclusion contains "B" (benign)
    vista_conclusion = line_object['VISta_conclusion']

    if 'B' in vista_conclusion:
        return True
    return False


def check_vista_conclusion_no_benign_or_vous(line_object):
    # Check if VISta conclusion does not contain "B" (benign) or "VOUS"
    vista_conclusion = line_object['VISta_conclusion']

    if 'B' not in vista_conclusion and 'VOUS' not in vista_conclusion:
        return True
    return False


def check_clinvar_conclusion_pathogenic(line_object):
    # Check if ClinVar conclusion contains "P" (pathogenic)
    clinvar_conclusion = line_object['ClinVar_conclusion']

    if 'P' in clinvar_conclusion:
        return True
    return False


def check_hgmd_conclusion(line_object):
    # Check if HGMD conclusion contains "DM" (disease-causing)
    hgmd_conclusion = line_object['HGMD_conclusion']

    if 'DM' in hgmd_conclusion:
        return True
    return False


def check_filter_and_zygosity(line_object):
    # Check filter value and zygosity value
    filter_value = line_object['filter']
    zygosity_value = line_object['zygosity']

    if filter_value == 'PASS' and zygosity_value != 'nocal':
        return True
    return False


def check_gnomad_maf_and_inheritance(line_object):
    # Check gnomAD highest subpopulation MAF and inheritance pattern
    inheritance_pattern = line_object['inheritance_pattern']
    gnomad_maf = float(line_object['gnomAD_highest_subpop_MAF'])

    if inheritance_pattern in ["AD", "AD/AR", "AR/AD", "XLR"]:
        if gnomad_maf < 0.5:
            return True
    elif inheritance_pattern == "AR" and gnomad_maf < 1:
        return True

    return False


def main(input_file, output_file):
    with open(input_file, 'r') as tsvfile, open(output_file, 'w', newline='') as tsv_output:
        tsv_reader = csv.DictReader(tsvfile, delimiter='\t')
        fieldnames = tsv_reader.fieldnames + ['Classification']

        tsv_writer = csv.DictWriter(tsv_output, fieldnames=fieldnames, delimiter='\t')
        tsv_writer.writeheader()

        for line in tsv_reader:
            filter_result = check_filter_and_zygosity(line)
            classification = 'Uncertain significance'  # Default classification

            if filter_result:
                gnomad_result = check_gnomad_maf_and_inheritance(line)
                hgmd_result = check_hgmd_conclusion(line)
                clinvar_result = check_clinvar_conclusion_pathogenic(line)
                vista_result = check_vista_conclusion_no_benign_or_vous(line)
                vista_curation_approved_result = check_vista_curation_status(line)
                vista_curation_needed_result = check_vista_curation_needed(line)
                vista_conclusion_benign_result = check_vista_conclusion_benign(line)

                if gnomad_result or hgmd_result or clinvar_result:
                    classification = 'Likely pathogenic'
                if vista_result and vista_curation_approved_result and not vista_curation_needed_result:
                    classification = 'Pathogenic'
                elif vista_result and not vista_curation_approved_result and not vista_curation_needed_result:
                    classification = 'Likely pathogenic'
                elif vista_conclusion_benign_result and vista_curation_approved_result:
                    classification = 'Benign'

            line['Classification'] = classification
            tsv_writer.writerow(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter and classify Single nucleotide variations (SNV).')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file')
    parser.add_argument('output_file', type=str, help='Path to the output TSV file')
    args = parser.parse_args()

    main(args.input_file, args.output_file)
