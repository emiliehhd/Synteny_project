import os
from Read_files import *
import subprocess

import psycopg2
conn = psycopg2.connect("dbname=base_hhed_projet")
cursor  = conn.cursor()


########## FONCTIONS INTERMEDIAIRES ##########
def cover(p_asembly, start, end):
    cursor.execute("SELECT longueur FROM genes WHERE id_gene = '%s'" % p_asembly)
    seq_len = cursor.fetchall()[0][0]
    return (end - start + 1) / seq_len

def verif_genomeNotInDataB(p_assembly):
    # retourne True si les genes du genomes ne sont pas dans la BDD
    not_done = True
    cursor.execute("SELECT DISTINCT assembly FROM genes")
    list_DL = [item[0] for item in cursor.fetchall()]

    if p_assembly in list_DL:
        not_done =False
    return not_done

########## FONCTIONS POUR L'INTERFACE GRAPHIQUE ##########
def guiCBBX_assembly():
    ## CREATION LISTE DE TOUS LES ASSEMBLY
    cursor.execute('SELECT assembly FROM genomes')
    rows_assembly   = cursor.fetchall()
    list_assembly   = [i[0] for i in rows_assembly]
    return list_assembly

def guiCBBX_blast():
    ## CREATION LISTE DES BLASTS DEJA REALISES
    cursor.execute('SELECT assembly_q, assembly_s FROM blast')
    rows_blast      = cursor.fetchall()
    list_blast      = [j[0] + " et " + j[1] for j in rows_blast]
    return  list_blast

def verif_BlastDone(entry1, entry2):
    """ Vérifie que le blast n'a pas été déjà fait (pour le messagebox error)
    Output : bool, indice du blast """
    blast_done = False
    res        = []

    cursor.execute("SELECT id_blast FROM blast "
                   "WHERE   (assembly_q = %s AND assembly_s = %s)"
                   "OR      (assembly_q = %s AND assembly_s = %s)",
                   (entry1, entry2, entry2, entry1))
    try:
        res = cursor.fetchall()
    except IOError:
        pass
    if res:
        blast_done = True
    return blast_done

def guiDynamicCCBX():
    regne_dict          = {}
    embranchement_dict  = {}
    classe_dict         = {}

    cursor.execute("SELECT DISTINCT regne FROM genomes")
    regne_fetch = cursor.fetchall()
    for r in regne_fetch:
        cursor.execute("SELECT DISTINCT embranchement FROM genomes "
                       "WHERE regne = '%s'" % r[0])
        regne_dict[r[0]] = [item[0]  for item in cursor.fetchall()]

    for i_embr in regne_dict.values():
        for e in i_embr:
            cursor.execute("SELECT DISTINCT classe FROM genomes "
                            "WHERE embranchement = '%s'" % e)
            embranchement_dict[e] = [item[0]  for item in cursor.fetchall()]

    for i_classe in embranchement_dict.values():
        for c in i_classe:
            cursor.execute("SELECT DISTINCT espece FROM genomes "
                           "WHERE classe = '%s'" % c)
            classe_dict[c] = [item[0] for item in cursor.fetchall()]

    return regne_dict, embranchement_dict, classe_dict

def recherche_taxo(p_regne, p_embr, p_classe, p_esp):
    requete = False
    list_assembly = []
    cursor.execute("SELECT souche, assembly FROM genomes "
                   "WHERE (regne = %s AND embranchement = %s AND classe = %s AND espece = %s)",
                   (p_regne, p_embr, p_classe, p_esp))
    try:
        list_assembly = ['souche :  '+item[0]+', assembly:  ' + item[1] for item in cursor.fetchall()]
    except IOError:
        pass
    if list_assembly:
        requete = True
    return requete, list_assembly

########## FONCTIONS POUR INSERT DANS DATABASE  ##########
def entreeDB_gene(faa_file):
    """ Entre les données des gènes d'un génome dans la table genes et le nombre de genes dans la table"""
    info_faa = read_faa(faa_file)

    for i in range(len(faa_file)):
        if faa_file[i:i+4] == "GCA_":
            assembly = faa_file[i : i+15]
            break

    for gene, value in info_faa.items():
        cursor.execute("INSERT INTO genes (id_gene, rang_gene, longueur, assembly) "
                       "VALUES ('%s', '%s', '%s', '%s')" %
                       (gene, value[0], value[1], assembly))

    conn.commit()


def entreeDB_blast_hits(blast_file):
    """ Input : nom du fichier blast out
        Entre les données des hits d'un blast dans la DB (modifie les tables blast et hits)"""
    info_blast  = read_blastOut(blast_file)
    last_hit    = ""
    assembly1, assembly2   = get_textAssembly(blast_file)

    ## Insertion dans la table blast : créé un identifiant pour le blast
    cursor.execute("INSERT INTO blast (id_blast, assembly_q, assembly_s)"
                   "VALUES (DEFAULT, '%s', '%s')" % (assembly1, assembly2))

    ## Récupération de l'identifiant du blast
    cursor.execute("SELECT id_blast FROM blast  "
                   "WHERE (assembly_q = '%s' AND assembly_s = '%s')" % (assembly1, assembly2))
    id_blast   = cursor.fetchall()[0][0]

    ## Insertions des hits du blast dans la table hits
    for hit in info_blast: #ordre dans hit : qseqid,  sseqid, pident, length, gapopen, qstart, qend, sstart, send, evalue
        Qcover = cover(hit[0], hit[5], hit[6])
        Scover = cover(hit[1], hit[7], hit[8])

        if last_hit != hit[0] + hit[1]:     # s'il y a deux hits ou plus pour le même query x subject, on garde le meilleur (le premier)
            cursor.execute("INSERT INTO hits (id_blast, query_id, subject_id, e_value, id_pourcent, start_subject, start_query, end_query, end_subject, cover_query, cover_subject) "
                           "VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')" % (id_blast, hit[0], hit[1], hit[9] , hit[2], hit[5], hit[6], hit[7],hit[8], Qcover, Scover))
            last_hit = hit[0] + hit[1]
    conn.commit()


########## FONCTIONS POUR RECUPERER DATA DANS DATABASE  ##########

def get_faaFile(entry, not_inDB):
    """ Récupère lien dans DB, télécharge et dézip si n'a jamais été fait, sinon retourne juste l'output
        output : nom du fichier faa """
    cursor.execute("SELECT genbank_ftp FROM genomes  "
                   "WHERE (assembly = '%s')"% entry)

    faa_file    = ''
    genbank     = []

    try:
        genbank = cursor.fetchall()
    except IOError:
        pass
    if genbank and not_inDB:  # cas où le fichier n'a jamais été téléchargé
        gb_link     = genbank[0][0]
        # link        = "http" + gb_link[3:] + "/" + gb_link[55:] + '_protein.faa.gz'
        link        = "http" + gb_link[3:] + "/" + gb_link[55:] + '_translated_cds.faa.gz'
        zip_file    = gb_link[55:] + '_translated_cds.faa.gz'
        faa_file    = gb_link[55:] + '_translated_cds.faa'
        ## Télécharment et unzip
        subprocess.run(["wget", link, "-O", zip_file])
        subprocess.run(["gunzip", zip_file])
        # subprocess.run(["rm "+ zip_file])

    if genbank and not not_inDB:
        gb_link = genbank[0][0]
        link = "http" + gb_link[3:] + "/" + gb_link[55:] + '_translated_cds.faa.gz'
        faa_file = gb_link[55:] + '_translated_cds.faa'
    return faa_file

def do_blast(faa_file1, faa_file2):
    output_blast = 'QUERY-{}__DB-{}.out'.format(faa_file1.strip('.faa'), faa_file2.strip('.faa'))

    subprocess.run(["blastp", "-query", faa_file1, "-subject", faa_file2, "-evalue", "1e-10", "-out", output_blast, "-outfmt",
                    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"])
    entreeDB_blast_hits(output_blast)
    # subprocess.run(["rm " + output_blast])


######### POUR LE DOTPLOT
def get_hits(entry, type_seuil, value_seuil):
    """ Retourne une liste contenant des listes [nom gene query, nom gene subject, rang gene query, rang gene subject]
    qui remplissent les conditions sélectionnées par l'utilisateur"""
    assembly1, assembly2 = get_textAssembly(entry)
    list_hits           = []                            # v1 avec id gene et position
    list_hitsPosition   = []                            # v2 avec jjuste les rangs
    len_genome1         = 0
    len_genome2         = 0
    id_blast            = 0
    col_seuil           = ""                            # nom de la colonne dans DB du seuil choisi
    operator            = ">"

    ## Récupération de la longueur des génomes : compte le nombre de gènes
    cursor.execute("SELECT COUNT(id_gene) FROM genes "
                   "WHERE assembly = '%s'"% assembly1)
    len_genome1 = cursor.fetchall()[0][0]
    cursor.execute("SELECT COUNT(id_gene) FROM genes "
                   "WHERE assembly = '%s'" % assembly2)
    len_genome2 = cursor.fetchall()[0][0]
    print("fin du count", len_genome1, len_genome2)
    cursor.execute("SELECT id_blast FROM blast "
                   "WHERE   (assembly_q = '%s' AND assembly_s = '%s')"
                   "OR      (assembly_q = '%s' AND assembly_s = '%s')"% (assembly1, assembly2, assembly2, assembly1))
    id_blast = cursor.fetchall()[0][0]

    if type_seuil == 0:                             # 0 : identite, 1 : evalue, 2, couverture
        col_seuil = "id_pourcent"
    elif type_seuil == 1:
        col_seuil = "e_value"
        operator  = "<"
    elif type_seuil == 2:
        col_seuil = "id_pourcent"

    print("debut select rang")

    print(col_seuil, operator, value_seuil, id_blast)
    cursor.execute("SELECT c.query_id, c.subject_id, a.rang_gene, b.rang_gene FROM hits as c "
                   "RIGHT JOIN genes AS a ON query_id = a.id_gene "
                   "RIGHT JOIN genes AS b ON subject_id = b.id_gene "
                   "WHERE ({} {} {} AND id_blast = {})".format(col_seuil, operator, value_seuil, id_blast))

    list_hits = [[i[0], i[1], i[2], i[3]] for i  in cursor.fetchall()]
    print(len(list_hits))

    print("fin de recuperation des rangs")

    return list_hits, len_genome1, len_genome2

## informations pour le dotplot
def get_infoGenomes(assembly):
    ## Recupere le rang taxonomique d'un génome
    cursor.execute("SELECT regne, embranchement, classe, espece, souche FROM genomes "
                   "WHERE assembly = '%s'" % assembly)
    taxonomie = cursor.fetchall()
    taxo_str    = "{}\n" \
                  "Règne                      : {}\n" \
                  "Embranchement    : {}\n" \
                  "Classe                      : {}\n" \
                  "Espèce                     : {}\n" \
                  "Souche                    : {}".format(assembly, taxonomie[0][0],
                                                 taxonomie[0][1], taxonomie[0][2], taxonomie[0][3], taxonomie[0][4])
    return taxo_str


