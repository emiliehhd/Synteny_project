import sys


######## READ FILES : BLAST ET FAA ########
def read_blastOut(nom_fi):
    # lecture d'un fichier blast
    # ordre out : qseqid,  sseqid, pident, length, gapopen, qstart, qend, sstart, send, evalue
    data_out = []
    try:
        file = open(nom_fi, "r")
    except IOError:
        print("Impossible d'ouvrir le fichier")
        sys.exit(1)

    for line in file:
        if line[0] != '#':
            data_l = line.split()
            data_query = data_l[0]
            data_subject = data_l[1]
            if data_query[:4] == "lcl|":       # pour les genes query avec "lcl|"
                data_query = data_query[4:]
            if data_subject[:4] == "lcl|":       # pour les genes subject avec "lcl|"
                data_subject = data_subject[4:]
            data_out.append([data_query, data_subject, float(data_l[2]), int(data_l[3]), int(data_l[5]), int(data_l[6]),
                             int(data_l[7]), int(data_l[8]), int(data_l[9]), float(data_l[10])])
    file.close()
    return data_out


def read_faa(nom_fi):
    # lecture d'un fichier fasta
    # output : dictionnaire avec (nom_gene : [rang, longueur])
    gene_info   = {}
    gene_name   = ''
    rang        = 1
    line_split  = ''
    try:
        file = open(nom_fi, "r")
    except IOError:
        print("Impossible d'ouvrir le fichier")
        sys.exit(1)

    for line in file:
        if line[0] == '>' and line[:5] != ">lcl|":
            gene_info[line.split()[0][1:]] = [rang, 0]
            gene_name = line.split()[0][1:]
            rang += 1
        if line[:5] == ">lcl|":
            line_split = line.split()
            gene_info[line_split[0][5:]] = [rang, 0]
            gene_name = line_split[0][5:]
            rang += 1
        if line[0] != '>':
            gene_info[gene_name][1] += len(line.strip())

    file.close()
    return gene_info


def get_textAssembly(string):
    """ Récupère les assembly dans une chaine de caractère"""
    assembly1 = ""
    assembly2 = ""
    for i in range(len(string)):
        if string[i:i+4] == "GCA_":
            if assembly1 == "":
                assembly1 = string[i : i+15]
            else:
                assembly2 = string[i : i+15]
    return assembly1, assembly2
