import sys
import tkinter as tk
from tkinter import messagebox as mb
from tkinter import ttk
# from DBconn import *
from DatabaseConn import *
from DotPlot import *
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


def seuil_OutOfRange(type_s, val_s):
    inRange = True
    val_s = float(val_s)
    if type_s == 0 and (val_s < 0 or val_s > 100):
        inRange = False
        mb.showinfo(title="Dotplot", message=" Le seuil d'identité n'est pas compris entre 0 et 100.")
    if type_s == 1 and (val_s < 0 or val_s > 1):
        inRange = False
        mb.showinfo(title="Dotplot", message=" Le seuil d'evalue n'est pas compris entre 0 et 1.")
    if type_s == 2 and (val_s < 0 or val_s > 100):
        inRange = False
        mb.showinfo(title="Dotplot", message=" Le seuil de couverture n'est pas compris entre 0 et 100.")
    return inRange


class Interface(tk.Tk):
    isInit = True
    regne, embranchement, classe = guiDynamicCCBX()

    def __init__(self):
        super().__init__()
        self.opt_var = None
        self.frame0             = None
        self.frame1             = None
        self.frame2             = None
        self.frame23            = None
        self.frame22            = None
        self.frame21            = None
        self.regne_ccbox        = None
        self.embranchement_ccbox = None
        self.classe_ccbox       = None
        self.espece_ccbx        = None
        self.stringence_scale   = None
        self.fenetre_scale      = None
        self.cover_spinbox      = None
        self.evalue_spinbox     = None
        self.ident_spinbox      = None
        self.var_critere        = None
        self.input1_cbbox       = None
        self.input2_cbbox       = None

        self.RootWindow_settings()
        self.Frames_settings()
        self.Frame_0()
        self.Frame_1()
        self.Frame_2()

        self.isInit = False

    def Lancer_recherche(self):
        var_regne   = self.regne_ccbox.get()
        var_embr    = self.embranchement_ccbox.get()
        var_classe  = self.classe_ccbox.get()
        var_esp     = self.espece_ccbx.get()

        corresp, list_correspAssembly = recherche_taxo(var_regne, var_embr, var_classe, var_esp)
        if not corresp:
            mb.showinfo(title="blast", message="Pas de genomes correspondant")
        else:
            self.Popup_taxo(list_correspAssembly, var_regne, var_embr, var_classe, var_esp)

    def Popup_taxo(self, l_correspAssembly, p_regne, pembr, p_classe, p_esp):
        popup = tk.Toplevel(self)
        popup.geometry("650x200")
        popup.title("Dotplot")
        popup['background'] = '#f2ede3'
        label_0 = ttk.Label(popup, text="Regne : {}, Embranchement : {}, Classe : {}, Espece : {}\n\n".format(p_regne, pembr, p_classe, p_esp), style="IG.TLabel")
        label_0.grid(row=0, column=0, columnspan=5, pady=(5, 0), padx=10)
        i=1
        for i_assembly in l_correspAssembly:
            b = ttk.Label(popup, text=i_assembly, style="IG.TLabel")
            b.grid(row=i, column=0, pady=(5, 0))
            i += 1

    def Lancer_blast(self):
        input1 = self.input1_cbbox.get()
        input2 = self.input2_cbbox.get()
        list_assembly = guiCBBX_assembly()

        ## Error : lance le blast sans entrées dans combobox query ou subject
        if (not input1 or not input2) and not self.isInit:
            mb.showinfo(title="blast", message="Selectionnez l'assembly de deux génomes")

        if input1 != '' and input2 != '':
            if verif_BlastDone(input1, input2):
                mb.showinfo(title="blast", message="Le blast pour ces deux génomes a déjà été réalisé !")
            elif input1 not in list_assembly and input2 not in list_assembly:
                mb.showinfo(title="blast", message="Les assembly saisis n'existent pas dans la base de données.")
            elif input1 not in list_assembly:
                mb.showinfo(title="blast", message="L'assembly du query n'existe pas dans la base de données.")
            elif input2 not in list_assembly:
                mb.showinfo(title="blast", message="L'assembly du subject n'existe pas dans la base de données.")
            else:
                input1_notInDB = verif_genomeNotInDataB(input1)
                input2_notInDB = verif_genomeNotInDataB(input2)
                input1_ffaFile = get_faaFile(input1, input1_notInDB)
                input2_ffaFile = get_faaFile(input2, input2_notInDB)
                print(input1_ffaFile)
                print(input2_ffaFile)

                ## Entrée des données
                if input1_notInDB:
                    entreeDB_gene(input1_ffaFile)
                if input2_notInDB:
                    entreeDB_gene(input2_ffaFile)

                print("\nLancement blast")
                do_blast(input1_ffaFile, input2_ffaFile)
                mb.showinfo(title="blast", message="Blast terminé")
                self.Frame_2()

    def Lancer_dotplot(self):
        var_NomBlast    = self.opt_var.get()
        type_seuil      = self.var_critere.get()        #0 : identite, 1 : evalue, 2, couverture
        val_seuil       = self.ident_spinbox.get()
        val_fenetre     = self.fenetre_scale.get()
        val_stringence  = self.stringence_scale.get()
        goBlast         = False

        if type_seuil == 1:
            val_seuil = self.evalue_spinbox.get()
        elif type_seuil == 2:
            val_seuil = self.cover_spinbox.get()

        ## ERROR : spinbox
        if val_seuil == '' and not self.isInit:
            list_typeSeuil = ["d'identité", "d'evalue", "de couverture"]
            mb.showinfo(title="Dotplot", message="Choisissez un seuil {}!".format(list_typeSeuil[type_seuil]))
        if val_seuil != "" and not self.isInit:
            goBlast = seuil_OutOfRange(type_seuil, val_seuil)

        if goBlast:
            list_hit1, list_hit2, len_genome1, len_genome2 = dotplot_final(var_NomBlast, type_seuil, val_seuil, val_fenetre, val_stringence)
            self.Popup_dotplot(var_NomBlast, list_hit1, list_hit2, len_genome1, len_genome2)

    def Popup_dotplot(self, blast_name, hits_query, hits_subject, len_query, len_subject):
        ## POP UP
        pop = tk.Toplevel(self)
        pop.geometry("1000x500")
        pop.title("Dotplot")
        pop['background'] = '#f2ede3'
        query, subject = get_textAssembly(blast_name)

        ## DOTPLOT AND CANVAS
        fig_dp, ax_dp = plt.subplots()
        ax_dp.plot(hits_subject, hits_query, marker=',', linestyle='None', color='black')
        ax_dp.set_title('Dotplot de {} et {}\n'.format(query, subject), fontsize=10, fontweight='bold',style='italic')
        ax_dp.set_ylim(0, len_query)
        ax_dp.set_xlim(0, len_subject)
        ax_dp.set_ylabel('query : ' + query + '\n',style='italic')
        ax_dp.set_xlabel('subject : ' + subject,style='italic')
        ax_dp.set_aspect('equal', adjustable="box")
        canvas = FigureCanvasTkAgg(fig_dp, master=pop)
        canvas.draw()

        ## INFOS GENOMES
        taxo_query      = get_infoGenomes(query)
        taxo_subject    = get_infoGenomes(subject)
        label_taxo      = ttk.Label(pop, text="Rangs taxonomiques des génomes ", style="IG.TLabel", font='Helvetica 10 bold')
        label_infos     = ttk.Label(pop, text=taxo_query+'\n\n\n'+taxo_subject, style="IG.TLabel")

        canvas.get_tk_widget().grid(row=0, column=0, rowspan=6, padx=5, pady=5,ipady=0, sticky="NSEW")
        label_taxo.grid(row=0, column=1)
        label_infos.grid(row=1, column=1, pady=(0, 0), padx=20)

    ## GUI STRUCTURE METHODS
    def RootWindow_settings(self):
        self['background'] = '#f2ede3'
        self.title("Projet 748")
        self.resizable()
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=10)

        ## STYLE
        style = ttk.Style()
        style.theme_use('default')
        style.configure('IG.TLabel', background="#f2ede3", foreground="#413936")
        style.configure('IG.TLabelframe', background="#f2ede3")
        style.configure('IG.TButton', foreground="#413936", background='#d3dbe1')
        style.configure('IG.TCombobox', foreground="#413936", selectbackground='#4682B4')
        style.configure('IG.TSpinbox', foreground='#413936', selectbackground='#4682B4')
        style.configure('IG.TRadiobutton', background='#f2ede3', foreground="#413936")
        style.configure('IG.TScale', background='#f2ede3', foreground="#413936", lightcolor="#413936")
        style.configure('IG.TMenubutton', background='#d3dbe1', foreground="#413936")

    def Frames_settings(self):
        self.frame0  = tk.LabelFrame(self, text = " Recherche ", background = '#f2ede3', foreground='#8B0000')
        self.frame1 = tk.LabelFrame(self, text = " 1) BLAST : sélection génomes ", background = '#f2ede3', foreground='#8B0000')
        self.frame2 = tk.LabelFrame(self, text = " 2) DOTPLOT ", background = '#f2ede3', foreground='#8B0000')
        self.frame21 = tk.LabelFrame(self.frame2, text = " Sélection BLAST ", background = '#f2ede3', foreground='#222e50')
        self.frame22 = tk.LabelFrame(self.frame2, text = " Conditions d'homologie", background = '#f2ede3', foreground='#222e50')
        self.frame23 = tk.LabelFrame(self.frame2, text = " Fenêtre et stringence", background = '#f2ede3', foreground='#222e50')

        self.frame0.grid(row=0, column=0, columnspan=2, padx=20, pady=10, ipady=7, sticky="NSEW")
        self.frame1.grid(row=1, column=0, columnspan=2, padx=20, pady=10, ipady=7, sticky="NSEW")
        self.frame2.grid(row=2, padx=20, pady=10, ipady=7, sticky="NSEW")
        self.frame21.grid(row=0, columnspan = 3, padx=10, pady=(10,5), ipady=7, sticky="NSEW")
        self.frame22.grid(row=1, columnspan = 3, padx=10, ipady=7, sticky="NSEW")
        self.frame23.grid(row=2, columnspan = 3, padx=10, pady=10, ipady=7, sticky="NSEW")

    def Frame_0(self):
        ttk.Label(self.frame0, text="Règne", style="IG.TLabel").grid(row=0, column=0)
        ttk.Label(self.frame0, text="Embranchement", style="IG.TLabel").grid(row=1, column=0)
        ttk.Label(self.frame0, text="Classe", style="IG.TLabel").grid(row=2, column=0)
        ttk.Label(self.frame0, text="Espèce", style="IG.TLabel").grid(row=3, column=0)


        self.espece_ccbx    = ttk.Combobox(self.frame0, width=40)
        self.classe_ccbox   = ttk.Combobox(self.frame0, width=40)
        self.classe_ccbox.bind('<<ComboboxSelected>>', self.Update_CCBX3)

        self.embranchement_ccbox = ttk.Combobox(self.frame0, width=40)
        self.embranchement_ccbox.bind('<<ComboboxSelected>>', self.Update_CCBX2)

        self.regne_ccbox = ttk.Combobox(self.frame0, width=40, values=list(self.regne.keys()))
        self.regne_ccbox.bind('<<ComboboxSelected>>', self.Update_CCBX1)

        button_taxo   = ttk.Button(self.frame0, text="Recherche", command=self.Lancer_recherche, style='IG.TButton', width=10)

        self.regne_ccbox.grid(row=0, column=1, padx=10)
        self.embranchement_ccbox.grid(row=1, column=1, padx=10)
        self.classe_ccbox.grid(row=2, column=1, padx=10)
        self.espece_ccbx.grid(row=3, column=1, padx=10)
        button_taxo.grid(row=2, column=2, padx=20)

    def Frame_1(self):
        ## LABEL
        input1_label        = ttk.Label(self.frame1, text="Query", style="IG.TLabel")
        input2_label        = ttk.Label(self.frame1, text="Subject", style="IG.TLabel")

        ## COMBOBOX
        self.input1_cbbox   = ttk.Combobox(self.frame1, values = guiCBBX_assembly(), style="IG.TCombobox")
        self.input2_cbbox   = ttk.Combobox(self.frame1, values = guiCBBX_assembly(), style="IG.TCombobox")
        button_blast        = ttk.Button(self.frame1, text="Lancer BLAST", command=self.Lancer_blast, style='IG.TButton', width=25)

        ## GRID
        input1_label.grid(row=0, column=0, pady=(5, 0))
        input2_label.grid(row=0, column=1)
        self.input1_cbbox.grid(row=1, column=0, padx=(15, 5))
        self.input2_cbbox.grid(row=1, column=1, padx=5)
        button_blast.grid(row=1, column=2,sticky="", padx=(10,5))

    def Frame_2(self):
        ## FRAME 2.1 : selection blast
        list_blastOptMenu   = guiCBBX_blast()
        BLAST_label         = ttk.Label(self.frame21, text="BLAST réalisés", style="IG.TLabel")
        self.opt_var = tk.StringVar()
        blast_optmenu  = ttk.OptionMenu(self.frame21,self.opt_var, list_blastOptMenu[0], *list_blastOptMenu, style="IG.TMenubutton")

        ## FRAME 2.2 : sélection des critères/conditions (threshold)
        criteres = ["Identité (%)", "Evalue", "Couverture"]
        self.var_critere = tk.IntVar()
        self.var_critere.set(0)
        i = 0
        for crit in criteres:
            c = ttk.Radiobutton(self.frame22, text=crit, value=i, variable=self.var_critere, style="IG.TRadiobutton")
            c.grid(row=0, column=i, pady=(5, 0))
            i += 1

        ident_label         = ttk.Label(self.frame22, text="Identité (0 à 100)", style="IG.TLabel")
        evalue_label        = ttk.Label(self.frame22, text="Evalue (0 à 1)", style="IG.TLabel")
        cover_label         = ttk.Label(self.frame22, text="Couverture (0 à 100)", style="IG.TLabel")
        self.ident_spinbox  = ttk.Spinbox(self.frame22, from_=0, to=100, increment=10, style='IG.TSpinbox')
        self.evalue_spinbox = ttk.Spinbox(self.frame22, format="%.5f", increment=0.00005, from_=0.000001, to=1, style='IG.TSpinbox')
        self.cover_spinbox  = ttk.Spinbox(self.frame22, from_=0, to=100, style='IG.TSpinbox')

        ## FRAME 2.3 : fenetre et stringence
        fenetre_label           = ttk.Label(self.frame23, text="Fenêtre : ", style="IG.TLabel")
        stringence_label       = ttk.Label(self.frame23, text="Stringence : ", style="IG.TLabel")
        self.fenetre_scale      = tk.Scale(self.frame23, from_=5, to=10, length=120, orient="horizontal", background = '#f2ede3')
        self.stringence_scale   = tk.Scale(self.frame23, from_=2, to=5, length=120, orient="horizontal", background = '#f2ede3')

        dotplot_button = ttk.Button(self.frame2, text="Afficher dotplot", command=self.Lancer_dotplot, style="IG.TButton")

        ## GRID
        BLAST_label.grid(row=0, column=0, pady=(10,0), padx=(10,7))
        blast_optmenu.grid(row=0, column=1, ipadx=40, pady=(10,0))
        ident_label.grid(row=1, column=0, pady=(5,0))
        self.ident_spinbox.grid(row=2, column=0, padx=(10,0))
        evalue_label.grid(row=1, column=1, pady=(10,0))
        self.evalue_spinbox.grid(row=2, column=1, padx=(15,0))
        cover_label.grid(row=1, column=2, pady=(10,0))
        self.cover_spinbox.grid(row=2, column=2, padx=(15,10))
        fenetre_label.grid(row=0, column=0, padx=(10,0), pady=(16,0))
        self.fenetre_scale.grid(row=0, column=1, padx=7)
        stringence_label.grid(row=0, column=2, padx=(20,7), pady=(16,0))
        self.stringence_scale.grid(row=0, column=3)
        dotplot_button.grid(row=3, column=0, columnspan=3, sticky="news", padx=20, pady=0)

    ## METHODES INTERMEDIARES

    def Update_CCBX1(self, event):
        self.embranchement_ccbox['values'] = self.regne[self.regne_ccbox.get()]

    def Update_CCBX2(self, event):
        self.classe_ccbox['values'] = self.embranchement[self.embranchement_ccbox.get()]

    def Update_CCBX3(self, event):
        self.espece_ccbx['values'] = self.classe[self.classe_ccbox.get()]

if __name__ == "__main__":
    interface = Interface()
    interface.mainloop()
