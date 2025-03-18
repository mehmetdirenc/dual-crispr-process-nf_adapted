import os
import sys


def read_main_library(zuber_lib_path):
    main_lib_dict = {}
    with open(zuber_lib_path, "r") as f:
        next(f)
        for line in f:
            entrez_id = line.split("\t")[4]
            if entrez_id + "_1" not in main_lib_dict:
                main_lib_dict[entrez_id + "_1"] = line
            else:
                main_lib_dict[entrez_id + "_2"] = line
    return main_lib_dict

def adapt_library(our_lib_path):
    lib_dict = {}
    with open(our_lib_path, "r") as f:
        next(f)
        for line in f:
            entrez_id = line.split("\t")[0]
            symbol = line.split("\t")[1]
            sequence = line.split("\t")[2]
            if sequence == "NA":
                continue
            sequence_matching = line.split("\t")[4]
            id = entrez_id + "_" + sequence + "_" + sequence_matching
            group = entrez_id
            context = "NA"
            gene_class = "NA"
            library_id = "dual_LM"
            adapted_line = (f"{id}\t{group}\t{sequence}\t{sequence_matching}\t"
                            f"{entrez_id}\t{symbol}\t{context}\t{gene_class}\t"
                            f"{library_id}\n")
            if entrez_id + "_1" not in lib_dict:
                lib_dict[entrez_id + "_1"] = adapted_line
            else:
                lib_dict[entrez_id + "_2"] = adapted_line
    return lib_dict


def combine_libraries(our_lib_path, main_lib_path, final_lib_path):
    our_lib_dict = adapt_library(our_lib_path)
    main_lib_dict = read_main_library(main_lib_path)
    with open(final_lib_path, "w") as f:
        f.write("id\tgroup\tsequence\tsequence_matching\tentrez_id\tsymbol\tcontext\tclass\tlibrary_id\n")
        for entrez_id in our_lib_dict:
            # if entrez_id not in main_lib_dict:
            f.write(our_lib_dict[entrez_id])
        for entrez_id in main_lib_dict:
            if main_lib_dict[entrez_id].replace("zuber_dual_mouse_library","dual_LM") in list(our_lib_dict.values()):
                continue
            if main_lib_dict[entrez_id] not in list(our_lib_dict.values()):
                f.write(main_lib_dict[entrez_id])
    return
if __name__ == '__main__':
    our_lib_path = "/home/direnc/inputs/kirsten/crispr/dual_LM_final.csv"
    main_lib_path = "/home/direnc/inputs/kirsten/crispr/Zuber_dual_mouse_library.txt"
    final_lib_path = "/home/direnc/inputs/kirsten/crispr/combined_library.txt"
    # final_lib_path = "/home/direnc/inputs/kirsten/crispr/dual_lm_adapted_library.txt"
    combine_libraries(our_lib_path, main_lib_path, final_lib_path)
    pass