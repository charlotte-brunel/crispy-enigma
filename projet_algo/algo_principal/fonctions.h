
/*******************************************************************************
 * * *                      FONCTIONS                                      * * *
 *******************************************************************************/

void importer_parametres(int* l, int* d, int* k, int* nb_masques);
void importer_sequences_fasta(TInfo_ensemble_sequences* ptr_info, TEnsemble_Sequences* ptr_ensemble );
void afficher_sequences(TInfo_ensemble_sequences* ptr_info, TEnsemble_Sequences* ptr_ensemble  );

int random_number(int max_number, int zero_excluded);
void generation_masque(int l, int* masque, int k);

void parcours_masque(int longueur_masque, int* adr_masque[longueur_masque], int nb_fenetre, int nb_sequence,
ptr_struct_seq* adr_tete_struct_sequence, TPtr_Cellkmer* adr_tete_liste_kmer,
TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos);

void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer,
TPtr_CellSequence* adr_liste_sequence, ptr_struct_seq* adr_liste_generation_sequence,
TPtr_CellPos* adr_liste_pos);


/*******************************************************************************
 * * *                   VARIABLES GLOBALES                                * * *
 *******************************************************************************/

extern bool convergence;
extern int l;
extern int d;
extern int k;
extern int nb_masques;
//
