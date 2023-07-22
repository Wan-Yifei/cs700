# Fliter rows for onetest HCNV 

# author: Yifei Wan (yifei.wan@sema4.com)

import sys
from var_check_module import *
from datetime import datetime


def get_kb_status(line_object, kb_db):
    """
    Query vonc_common_data.vonc_cnv_kb_current
    to get the latest value of KB classification of CNV
    """
    ## Query the database
    records = kb_db.kb_query(line_object.fileds["call_type"], line_object.transcript, line_object.fileds["exon_start"], line_object.fileds["exon_end"]) 
    if not records:
        return None
    else:
        kb_status = {}
        kb_status["classification_nm"] = records[0][0]
        kb_status["transaction_start_tm"] = records[0][1]
        kb_status["transaction_until_tm"] = records[0][2]
        return kb_status

def get_hc_cnv_status(line_object, kb_db):
    """
    Query the vonc_hcancer.var_hc_cnv_data
    to get the classification of CNV in other previous cases
    """
    ## Query the database
    records = kb_db.kb_query_hc(line_object.fileds["call_type"], line_object.transcript, line_object.fileds["exon_start"], line_object.fileds["exon_end"]) 
    if not records:
        return None
    else:
        hc_cnv_status = {}
        hc_cnv_status["patient_cnv_patho_clssfctn_nm"] = records[0][0]
        hc_cnv_status["transaction_start_tm"] = records[0][1] 
        hc_cnv_status["transaction_until_tm"] = records[0][2]
        return hc_cnv_status

def check_qc_filter(line_object):
    if line_object.filter == "PASS":
        return True
    else:
        return False
    
def check_unclassfied(line_object): ## filter 1
    if line_object.fileds["classification"] == "Unclassified" and check_qc_filter(line_object):
        return True
    else:
        return False

def check_kb_benign(kb_status): ## filter 2
    """
    The the latest classification of CNV is B/LB in KB. 
    """
    if not kb_status:
        return False
    current = datetime.utcnow()
    if kb_status["transaction_start_tm"] <= current and kb_status["transaction_until_tm"] > current \
        and "B" in kb_status["classification_nm"]:
        return True
    else:
        return False   
    
def check_hc_cnv_benign(hc_cnv_status): ## filter 3
    """
    Check the CNV in previous cases having approved (in KB) B/LB.
    """
    if not hc_cnv_status:
        return False
    current = datetime.utcnow()
    if hc_cnv_status["transaction_start_tm"] <= current and hc_cnv_status["transaction_until_tm"] > current \
        and "B" in hc_cnv_status["patient_cnv_patho_clssfctn_nm"]:
        return True
    else:
        return False   

def check_kb_patho(kb_status, line_object): ## filter 4
    """
    the latest classification of cnv is P/LP in KB.
    """
    if not kb_status:
        return False
    current = datetime.utcnow()
    if kb_status["transaction_start_tm"] <= current and kb_status["transaction_until_tm"] > current \
        and "P" in kb_status["classification_nm"] and check_qc_filter(line_object): 
        return True
    else:
        return False   

def check_kb_vus(kb_status, line_object): ## filter 5
    """
    the latest classification of cnv is VOUS in KB. Not expired.
    """
    if not kb_status:
        return False
    current = datetime.utcnow()
    if kb_status["transaction_start_tm"] <= current and kb_status["transaction_until_tm"] > current \
        and "V" in kb_status["classification_nm"] and \
         (current - kb_status["transaction_start_tm"]).days/30 < 12 and check_qc_filter(line_object):
        return True
    else:
        return False

def check_kb_expired_vus(kb_status, line_object): ## filter 6
    """
    The latest classification of CNV is VUS but expired.
    """
    if not kb_status:
        return False
    current = datetime.utcnow()
    if kb_status["transaction_start_tm"] <= current and kb_status["transaction_until_tm"] > current \
        and "V" in kb_status["classification_nm"] and (current - kb_status["transaction_start_tm"]).days/30 >= 12 \
        and check_qc_filter(line_object):
        return True
    else:
        return False   

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    config_file = "config_hcancer.json"
    get_credentials(config_file)
    kb_db = KB_record(config_file) ## engine to the database
    with open(input_file, "r") as raw, open(output_file, "w+") as output:
        line = raw.readline()
        header = line
        print >> output, header.strip()
        line = raw.readline()
        while line: 
            line_object = Variant_prepare(header, line, kb_db) ## the object of one line
            kb_status = get_kb_status(line_object, kb_db)
            hc_cnv_status = get_hc_cnv_status(line_object, kb_db)
            if check_unclassfied(line_object): ## filter 1. Unclassified (filter 3 in sql code)
                line_object.set_select(1) ## select the variant 
            if check_kb_benign(kb_status): ## filter 2. check B/LB in KB (filter 1 in sql code)
                line_object.set_select(0) ## not select 
            #if check_hc_cnv_benign(hc_cnv_status): ## filter 3. check B/LB in old cases (filter 5 in sql code)
            #    line_object.set_select(0) ## not select 
            if check_kb_patho(kb_status, line_object): ## filter 4. P/LP in KB (filter 2 in sql code)
                line_object.set_select(1) ## select the variant 
                line_object.set_report(1) ## report the variant
                line_object.patient_case_variant_curation_action_id = "11" ## When equals 11, appers variant in the report candidate table on UI
            if check_kb_vus(kb_status, line_object): ## filter 5. VUS in KB not expired (filter 2 in sql code)
                line_object.set_select(1) ## select the variant 
                line_object.set_report(1) ## report the variant
                line_object.patient_case_variant_curation_action_id = "11"
            if check_kb_expired_vus(kb_status, line_object): ## filter 6. VOUS in KB, expired (filter 4 in sql code) 
                line_object.set_select(1) ## select the variant 

            new_line = "\t".join(line.strip().split("\t")[:len(line.strip().split("\t")) - 4]) + "\t" + "%d"%line_object.is_selected_in + "\t" + "%s"%line_object.is_report_candidate_in + "\t"  + "%s"%line_object.actor_user_id + "\t" + "%s"%line_object.patient_case_variant_curation_action_id
            del line_object 

            #print(new_line, file = output)
            print >> output, new_line
            line = raw.readline() 
        kb_db.kb_close()
            
if __name__ == "__main__":
    main()                    
