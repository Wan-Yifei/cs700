# Variant class for variant filter
# author: Yifei.wan

import copy
import sys
from db_client import Database_visit
from datetime import datetime

class Variant(object):
    def __init__(self, header, line):
        cols = header.strip().split("\t")
        self.fileds = {col:line.strip().split("\t")[cols.index(col)] for col in cols} 
        self.filter = self.fileds["filter"]
        try:
            self.vkey = self.fileds["vkey"]
        except:
            pass
        try:
            self.gene = self.fileds["gene"]
        except:
            pass
        try:
            self.transcript = self.fileds["transcript"] 
        except:
            pass
        try:
            self.zygosity = self.fileds["zygosity"]
        except:
            pass
        try:
            self.hgmd = self.fileds["HGMD_conclusion"]
        except:
            pass
        try:
            self.clinvar = self.fileds["ClinVar_conclusion"] 
        except:
            pass
        try:
            self.inheritance_pattern = self.fileds["inheritance_pattern"]
        except:
            pass
        try:
            self.highest_subpop_MAF = self.fileds["gnomAD_highest_subpop_MAF"]
        except:
            pass
        try:
            self.effect = self.fileds["effect"]
        except:
            pass
        try:
            self.roi_candidate = int(self.fileds["ROI_candidate"])
        except:
            pass

class Variant_prepare(Variant):
    def __init__(self, header, line, vista_db):
        Variant.__init__(self, header, line)
        self.is_selected_in = 0
        self.is_report_candidate_in = 0
        self.actor_user_id = 4
        self.patient_case_variant_curation_action_id = 0
        self.export_eligibility_group_id = 1
        self.kb = None
        self.vista = None
        self.mul = None
        self.rule = [] ## what are filters passed?

    def set_export_eligibility_group_id(self, value):
        self.export_eligibility_group_id = value

    def set_select(self, value):
        self.is_selected_in = value

    def set_report(self, value):
        self.is_report_candidate_in = value

    def set_patient_case_curation_action(self, value):
        self.patient_case_variant_curation_action_id = value

    def hgmd_check(self):
        if "DM" in self.hgmd:
            self.hgmd_flag = True ## pathogenic
        else:
            self.hgmd_flag = False ## benign        

    def get_kb_status(self, kb_record):
        current = datetime.utcnow()
        #kb_flag = 0 ## init status, 0 = failed.
        records = kb_record.kb_query(self.fileds["call_type"], self.transcript, self.fileds["exon_start"], self.fileds["exon_end"]) 
        if not records:
            return None
        else:
            kb_status = {}
            kb["classification_nm"] = record[0][0]
            kb["transaction_start_tm"] = record[0][1]
            kb["transaction_until_tm"] = record[0][2]
            return kb_status

    def get_hc_cnv_status(self, kb_record):
        # Check data in var_hc_cnv_data
        current = datetime.utcnow()
        #kb_flag = 0 ## init status, 0 = failed.
        records = kb_record.kb_query_hc(self.fileds["call_type"], self.transcript, self.fileds["exon_start"], self.fileds["exon_end"]) 
        if not records:
            return None
        else:
            hc_cnv_status = {}
            hc_cnv_status["patient_cnv_patho_clssfctn_nm"] = records[0][0]
            hc_cnv_status["transaction_start_tm"] = records[0][1] 
            hc_cnv_status["transaction_until_tm"] = records[0][2]
            return hc_cnv_status

    def get_vista_status(self, vista_record):
        current = datetime.utcnow()
        vista_flag = 0 ## init status, 0 = failed.
        records = vista_record.vista_query(self.vkey, self.gene, self.transcript) 
        if not records:
            return None
        vista_conclusion_nm, vista_curation_status_nm, transaction_start_tm, transaction_until_tm = records[0] 
        if "B" in vista_conclusion_nm:
            vista_flag += 1
        if "VOUS" in vista_conclusion_nm:
            vista_flag += 10
        if "P" in vista_conclusion_nm:
            vista_flag += 100
        if transaction_start_tm <= current and  transaction_until_tm > current:
            vista_flag += 1000
        if (current - transaction_start_tm).days/30 > 12:
            vista_flag += 10000
        if vista_curation_status_nm == "APPROVED":
            vista_flag += 100000 ## Approved
        if vista_curation_status_nm == "Precalculated":
            vista_flag += 1000000 ## Not approved
        self.vista = int("%d"%vista_flag, 2)

    def mul_cnv_status(self, mul_cnv_record):
        current = datetime.utcnow()
        mul_cnv_flag = 0
        record = mul_cnv_record.mul_cnv_query_hc(self.fileds["call_type"], self.fileds["chrom"], self.fileds["segL"], self.fileds["segR"], self.fileds["segLen"]) 
        if not record:
            self.mul = int("%d"%mul_cnv_flag, 2) ## Not in bk 
        else:
            curation_status_nm, patient_cnv_patho_clssfctn_nm = record[0] 
            current = datetime.utcnow()
            if curation_status_nm == "in KB":
                mul_cnv_flag += 1
            if "B" in patient_cnv_patho_clssfctn_nm:
                mul_cnv_flag += 10
        
        ## check UI table
        record = mul_cnv_record.mul_cnv_query_ui(self.fileds["call_type"], self.fileds["chrom"], self.fileds["segL"], self.fileds["segR"], self.fileds["segLen"]) 
        if not record:
            self.mul = int("%d"%mul_cnv_flag, 2) ## Not in bk 
            return
        curation_status_nm, patient_cnv_patho_clssfctn_nm, kb_curation_tm = record[0]
        if curation_status_nm == "in KB":
            mul_cnv_flag += 100 
        if "V" in patient_cnv_patho_clssfctn_nm:
            mul_cnv_flag += 1000 ## UI class
        if (current - kb_curation_tm).days/30 > 12:
            mul_cnv_flag += 10000
        self.mul = int("%d"%mul_cnv_flag, 2)
                    
    def freq_check(self):

        if self.inheritance_pattern in ["AD", "AD/AR", "AR/AD", "XLR"]:
            if float(self.highest_subpop_MAF) < 0.5:
                self.freq_flag = True ## pass
            else:
                self.freq_flag = False ## failed
        elif self.inheritance_pattern == "AR":
            if float(self.highest_subpop_MAF) < 1:
                self.freq_flag = True ## pass
            else:
                self.freq_flag = False ## failed
        else: 
            self.freq_flag = False ## No recognized pattern

    def effect_check(self):
        accept_effects = ["frameshit", "stopgain", "stop_gain", "stop_lost", "stoploss", "start_lost", "startloss", "splice_donor", \
                            "splice donor", "splice_acceptor", "splice acceptor", "initiator_codon_variant", "initiator codon", "missense" \
                            "synonymous"]
        for eff in accept_effects:
            if eff in self.effect:
                self.effect_flag = True ## pass
                return
            else:
                self.effect_flag = False ## not a required effect

    def region_check(self):
        """
        gene based ROI special rules.
        """
        if self.roi_candidate == 1:
            self.roi_candidate_flag = True ## in intersted region
        else:
            self.roi_candidate_flag = False


class Vista_record(object): 
    def __init__(self, config):
        self.db = Database_visit(config)
        self.db.connect() ## connect to database
        self.db.open_cursor() ## open a cursor
        
    def vista_query(self, vkey, gene, transcript):
        command = """select VISta_classification, status, mvl_last_updated_dt, created_dt from physical_view_var_pending\
         where vkey = "%s" and symbol = "%s" and substring_index(refseq, ".", 1) = substring_index("%s", ".", 1);"""%(vkey, gene, transcript)
        self.db.query(command)
        if self.db.cursor.rowcount == 0:
            return None
        else:
            return self.db.cursor.fetchall()

    def vista_close(self):
        self.db.close_cursor()
        self.db.close_database()

class KB_record(object): 
    def __init__(self, config):
        self.db = Database_visit(config)
        self.db.connect() ## connect to database
        self.db.open_cursor() ## open a cursor
        
    def kb_query(self, call_type, transcript, exon_start, exon_end):
        command = """select classification_nm, transaction_start_tm, transaction_until_tm \
from vonc_common_data.vonc_cnv_kb_current \
where cnv_call_type_nm = "%s" and transcript_nm = "%s" and exon_start_num = %s and \
exon_end_num = %s;"""%(call_type, transcript, exon_start, exon_end)
        self.db.query(command)
        if self.db.cursor.rowcount == 0:
            return None
        else:
            return self.db.cursor.fetchall()

    def kb_query_hc(self, call_type, transcript, exon_start, exon_end):
        command = """select patient_cnv_patho_clssfctn_nm, transaction_start_tm, transaction_until_tm \
from var_hc_cnv_data \
where cnv_call_type_nm = "%s" and transcript_nm = "%s" and exon_start_num = %s and \
exon_end_num = %s and curation_status_nm = "in KB";"""%(call_type, transcript, exon_start, exon_end)
        self.db.query(command)
        if self.db.cursor.rowcount == 0:
            return None
        else:
            return self.db.cursor.fetchall()

    def kb_close(self):
        self.db.close_cursor()
        self.db.close_database()

class Multi_cnv_record(object):
    def __init__(self, config):
        self.db = Database_visit(config)
        self.db.connect() ## connect to database
        self.db.open_cursor() ## open a cursor

    def mul_cnv_query_hc(self, call_type, chrom, segL, segR, segLen):
        command = """select patient_cnv_patho_clssfctn_nm, kb_curation_tm \
from var_hc_cnv_multi_gene_data \
where cnv_call_type_nm = "%s" and seg_chrom_nm = "%s" and seg_start_pos_num = "%s" and \
seg_end_pos_num = "%s" and seg_length_am = "%s" and curation_status_nm = "in KB";"""%(call_type, chrom, segL, segR, segLen)
        self.db.query(command)
        if self.db.cursor.rowcount == 0:
            return None
        else:
            return self.db.cursor.fetchall()

    def mul_cnv_query_ui(self, call_type, chrom, segL, segR, segLen):
        command = """select patient_cnv_patho_clssfctn_nm, kb_curation_tm \
from v_hc_cnv_ui \
where cnv_call_type_nm = "%s" and seg_chrom_nm = "%s" and seg_start_pos_num = "%s" and \
seg_end_pos_num = "%s" and seg_length_am = "%s" and curation_status_nm = "in KB";"""%(call_type, chrom, segL, segR, segLen)
        self.db.query(command)
        if self.db.cursor.rowcount == 0:
            return None
        else:
            return self.db.cursor.fetchall()

    def mul_cnv_vus_query_kb(self, call_type, chrom, segL, segR):
        command = """select classification_nm, transaction_start_tm from vonc_common_data.vonc_multi_gene_cnv_kb_current where cnv_call_type_nm = "%s" and chrom_nm = "%s" and cnv_genomic_pos_start_num = "%s" and cnv_genomic_pos_end_num = "%s";"""%(call_type, chrom, segL, segR)
        self.db.query(command)
        if self.db.cursor.rowcount == 0:
            return None
        else:
            return self.db.cursor.fetchall()

    def mul_close(self):
        self.db.close_cursor()
        self.db.close_database()
