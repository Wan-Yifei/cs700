# Fliter rows for onetest SV

# author: Yifei Wan

import sys
import json
from var_check_module import * 
from datetime import datetime

def qc_check(line_dict):
    if line_dict["filter"] == "PASS":
        return True
    else:
        return False

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    with open(input_file, "r") as raw, open(output_file, "w+") as output:
        line = raw.readline()
        header = line
        cols = header.strip().split("\t")
        #print(header.strip(), file = output)
        print >> output, header.strip()
        line = raw.readline()
        while line: 
            fileds = {col:line.strip().split("\t")[cols.index(col)] for col in cols}
            #var.mapping_quality = 30
            #if fileds["filter"] == "PASS":
            #    flag = True
            #else:
            #    flag = False
            flag = qc_check(fileds) 
            if flag:
                update_flag = line.strip().split("\t")
                update_flag[-4:] = ["1", "1", "4", "11"]
                new_line = "\t".join(update_flag) 
                line = new_line.strip()
            #print(new_line, file = output)
            line = line.strip()
            print >> output, line
            line = raw.readline() 


if __name__ == "__main__":
    main()
