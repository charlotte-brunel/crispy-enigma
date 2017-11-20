/*
 *******************************
 * STRUCTURE DICTIONNAIRE KMER *
 *******************************
1er élément: liste chainée des kmers

*/
typedef struct TCellkmer TCellkmer;
struct TCellkmer{
    struct TCellkmer* suiv_kmer; //Pointer sur le kmer suivant
    struct TCellSequence* tete_sequence; //Pointer qui pointe sur une liste chainée de séquence où le kmer est retrouvé
    char kmer[]; //"TCC" par exemple

};
typedef TCellkmer* TPtr_Cellkmer; //Pointeur sur Tcellkmer

//deuxième élément du dictionnaire kmer: liste chainée de sequence possédant le kmer:

typedef struct TCellSequence TCellSequence;
struct TCellSequence{
    int sequence; //sequence 1 --> 1; sequence 2 --> 2
    struct TCellSequence* suiv_sequence; // Pointeur qui pointe l'élément suivant de la liste chainée de séquence
    struct TCellPos* tete_pos; // Pointeur qui pointe sur le premier élément de la structure position qui répertorie toutes les positions où le kmer a été trouvé dans une séquence
};
typedef TCellSequence* TPtr_CellSequence; // Pointeur sur TCellSéquence

// troisième élément du dictionnaire de kmer: liste chainée de position où le kmer a été trouvé dans une séquence.

typedef struct TCellPos TCellPos;
struct TCellPos{
    int position; // position du kmer dans la séquence
    struct TCellPos* suiv_pos; //pointe sur la position suivante
};
typedef TCellPos* TPtr_CellPos;

//NOormalement ça tu n'en as pas besoin c'est la structure qui contient les séquences qu"on étudie
/* ***********************************
   * STRUCTURE CHAINEE DES SEQUENCES *
   ***********************************


typedef struct structure_sequence structure_sequence;
struct structure_sequence{
	int numero_sequence; // 1 (indication de la séquence)
	char sequence[31]; // 'ATCGGACG...' séquencede nucléotide
	struct structure_sequence *next_sequence; //pointeur sur le prochain numéro de séquence (ex: '>Seq2:')
	};

typedef structure_sequence* ptr_struct_seq; // on créé le type ptr_struct_seq qui est un pointeur sur la structure de séquence.
*/




void parcours_masque(int longueur_masque, int* adr_masque[longueur_masque], int nb_fenetre, int nb_sequence, ptr_struct_seq* adr_tete_struct_sequence, TPtr_Cellkmer* adr_tete_liste_kmer, TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos){
    ptr_struct_seq p_generation_seq = *adr_tete_struct_sequence;
    TPtr_Cellkmer p_kmer= *adr_tete_liste_kmer;
    TPtr_CellSequence p_sequence= *adr_tete_liste_sequence;
    TPtr_CellPos p_pos= *adr_tete_liste_pos;
    int position=0;
    int cpt_pos_masque, position_kmer;
    int pos_kmer =0;
    char k_mer[nb_fenetre+1]; // le k_mer mesure la taille du nombre de fenêtre ouverte dans le masque
    while (p_generation_seq != NULL){
        while (position<30){
                for(cpt_pos_masque=0; cpt_pos_masque< longueur_masque; cpt_pos_masque++){ // on parcourt les nucleotides sous le masque
                    if (adr_masque[cpt_pos_masque]==1){ //si la fenêtre du masque est ouverte
                        k_mer[pos_kmer]= p_generation_seq->sequence[position];

                        if (pos_kmer==0){
                            position_kmer=position;
                        }
                        pos_kmer++;
                        if (pos_kmer==nb_fenetre){
                            k_mer[pos_kmer]= '\0' ;
                            pos_kmer=0;
                        }

                    }
                    position= position +1;
                }
            generation_kmer(position_kmer, k_mer, &p_kmer, &p_sequence, &p_generation_seq, &p_pos);
        }
        position=0;
        p_generation_seq= p_generation_seq->next_sequence;
    }
    return;
}


void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer, TPtr_CellSequence* adr_liste_sequence, ptr_struct_seq* adr_liste_generation_sequence, TPtr_CellPos* adr_liste_pos){
TPtr_Cellkmer temp_p_kmer= *adr_liste_kmer; //on créé des pointeurs temporaire pour parcourir les listes
ptr_struct_seq p_generation_seq= *adr_liste_generation_sequence;
if (temp_p_kmer->suiv_kmer==NULL){ //cas du premier element de la liste
    TPtr_Cellkmer nouveau_kmer=(TCellkmer*)malloc(sizeof(TCellkmer));
    strcpy(temp_p_kmer->kmer, k_mer);
    nouveau_kmer->suiv_kmer= NULL;
    TPtr_CellSequence nouvelle_sequence=(TCellSequence*)malloc(sizeof(TCellSequence));
    temp_p_kmer->suiv_kmer= nouveau_kmer;
    temp_p_kmer->tete_sequence=nouvelle_sequence;
    nouvelle_sequence->sequence=p_generation_seq->numero_sequence;
    nouvelle_sequence->suiv_sequence= NULL;
    TPtr_CellPos nouvelle_position=(TCellPos*)malloc(sizeof(TCellPos));
    nouvelle_sequence->tete_pos= nouvelle_position;
    nouvelle_position->position=position_kmer;
    nouvelle_position->suiv_pos= NULL;
    return;
    }else{
        while (temp_p_kmer->suiv_kmer != NULL){
            if (strcmp(temp_p_kmer->kmer, k_mer)==0){ //si le kmer a déjà été trouvé
                TPtr_CellSequence p_liste_sequence=temp_p_kmer->tete_sequence;
                //premier element de la liste chainee de sequence:
                while (p_liste_sequence->suiv_sequence != NULL){
                    if (p_generation_seq->numero_sequence==p_liste_sequence->sequence){ //Si le kmer a déjà été trouvé dans cette séquence
                        TPtr_CellPos p_liste_pos= p_liste_sequence->tete_pos;
                        while (p_liste_pos->suiv_pos != NULL){
                            p_liste_pos=p_liste_pos->suiv_pos;
                            }
                        TPtr_CellPos nouvelle_position= (TCellPos*)malloc(sizeof(TCellPos)); // On créé une nouvelle brique de pôsition
                        nouvelle_position->position=position_kmer;
                        nouvelle_position->suiv_pos=NULL;
                        p_liste_pos->suiv_pos= nouvelle_position;
                        return;
                        }
                    p_liste_sequence=p_liste_sequence->suiv_sequence;
                    }
                if (p_generation_seq->numero_sequence==p_liste_sequence->sequence){ //Si le kmer a déjà été trouvé dans cette séquence
                    TPtr_CellPos p_liste_pos= p_liste_sequence->tete_pos;
                    while (p_liste_pos->suiv_pos != NULL){
                            p_liste_pos=p_liste_pos->suiv_pos;
                            }
                    TPtr_CellPos nouvelle_position= (TCellPos*)malloc(sizeof(TCellPos)); // On créé une nouvelle brique de pôsition
                    nouvelle_position->position=position_kmer;
                    nouvelle_position->suiv_pos=NULL;
                    p_liste_pos->suiv_pos= nouvelle_position;
                    return;
                    }
                TPtr_CellSequence nouvelle_sequence= (TCellSequence*)malloc(sizeof(TCellSequence));
                nouvelle_sequence->sequence= p_generation_seq->numero_sequence;
                nouvelle_sequence->suiv_sequence=NULL;
                p_liste_sequence->suiv_sequence= nouvelle_sequence;
                TPtr_CellPos nouvelle_tete_pos= (TCellPos*)malloc(sizeof(TCellPos)); // on créé une nouvelle tete de liste de position
                nouvelle_sequence->tete_pos= nouvelle_tete_pos;
                nouvelle_tete_pos->position=position_kmer;
                nouvelle_tete_pos->suiv_pos= NULL;
                return;

            }
            temp_p_kmer=temp_p_kmer->suiv_kmer;
        }



       // Si ce kmer n'a pas encore été trouvé: ajout en fin de boucle
        TPtr_Cellkmer nouveau_kmer=(TCellkmer*)malloc(sizeof(TCellkmer));
        strcpy(temp_p_kmer->kmer, k_mer);
        nouveau_kmer->suiv_kmer= NULL;
        TPtr_CellSequence nouvelle_sequence=(TCellSequence*)malloc(sizeof(TCellSequence));
        temp_p_kmer->suiv_kmer= nouveau_kmer;
        temp_p_kmer->tete_sequence=nouvelle_sequence;
        nouvelle_sequence->sequence=p_generation_seq->numero_sequence;
        nouvelle_sequence->suiv_sequence= NULL;
        TPtr_CellPos nouvelle_position=(TCellPos*)malloc(sizeof(TCellPos));
        nouvelle_sequence->tete_pos= nouvelle_position;
        nouvelle_position->position=position_kmer;
        nouvelle_position->suiv_pos= NULL;
        return;
    }

}





//APPEL DE LA FONCTION:

//Variables:
//TPtr_Cellkmer element_kmer=malloc(sizeof(TCellkmer));
//TPtr_Cellkmer tete_liste_kmer= element_kmer;

//TPtr_CellSequence element_sequence= malloc(sizeof(TCellSequence));
//TPtr_CellSequence tete_liste_sequence= element_sequence;

//TPtr_CellPos element_pos= malloc(sizeof(TCellPos));
//TPtr_CellPos tete_liste_pos= element_pos;

//ptr_struct_seq element_generation_sequence=malloc(sizeof(structure_sequence)); // On créé le premier élément de la structure
//ptr_struct_seq tete_liste_pour_parcours_masque= element_generation_sequence;

//parcours_masque(longueur_masque, &adr_masque[longueur_masque], nb_fenetre, nb_sequence, &tete_liste_pour_parcours_masque, &tete_liste_kmer, &tete_liste_sequence, &tete_liste_pos);
// &tete_liste_pour parcours_masque c'est un pointeur sur une liste chainée de séquence (qu'on récupère dans le dico)
