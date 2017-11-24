/* A FAIRE
-r�tuiliser la fonction insert motif pour tester 1 premier calcul de PSSM
-faire en sorte que la generation de masque est toujours au minimum deux fenetres visible et au maximum quatre.

*/

#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<string.h>
#include "function.h"


int main()
{

FILE *file_info=NULL; //fichier contenant les infos sur les s�quences
int longueur_masque = 5;
int masque[longueur_masque];
memset(masque, 0, sizeof masque);
//int* p_masque= masque;
int nb_fenetre=2;
TPtr_Cellkmer element_kmer=malloc(sizeof(TCellkmer));
TPtr_Cellkmer tete_liste_kmer= element_kmer;
TPtr_Cellkmer tete_liste_kmer2= element_kmer;
TPtr_Cellkmer tete_liste_kmer3= element_kmer;
TPtr_Cellkmer tete_liste_kmer4= element_kmer;
element_kmer->suiv_kmer=NULL;
TPtr_CellSequence element_sequence= malloc(sizeof(TCellSequence));
TPtr_CellSequence tete_liste_sequence= element_sequence;
TPtr_CellSequence tete_liste_sequence2= element_sequence;
TPtr_CellSequence tete_liste_sequence3= element_sequence;
element_sequence->suiv_sequence=NULL;
TPtr_CellPos element_pos= malloc(sizeof(TCellPos));
TPtr_CellPos tete_liste_pos= element_pos;
TPtr_CellPos tete_liste_pos2= element_pos;
TPtr_CellPos tete_liste_pos3= element_pos;
element_pos->suiv_pos=NULL;
TPtr_Cellkmer_selectionne element_liste_kmer_selectionne= malloc(sizeof(TCellkmer_selectionne));
TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne= element_liste_kmer_selectionne;
TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne2= element_liste_kmer_selectionne;
TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne_pour_calcul= element_liste_kmer_selectionne;
TPtr_Cell_Motif_PSSM element_liste_motif_PSSM= malloc(sizeof(TCell_Motif_PSSM));
TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM= element_liste_motif_PSSM;
TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM2= element_liste_motif_PSSM;
TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM_pour_calcul= element_liste_motif_PSSM;
element_liste_motif_PSSM->suiv_motif=NULL;

//MATRICE PSSM ET MOTIF CONSENSUS:

double **ptr_matrice_PSSM;


/* GENERATION ALEATOIRE DE DIX SEQUENCES POUR TEST */
//D�claration des variables
srand(time(0));

int nb_sequence= 10; // nombre de s�quence � g�n�rer
int i,j; //it�rateur de boucle
char nucleo= ' ';

ptr_struct_seq element_generation_sequence=malloc(sizeof(structure_sequence)); // On cr�� le premier �l�ment de la structure
ptr_struct_seq tete_liste_pour_parcours_masque= element_generation_sequence;
ptr_struct_seq tete_liste_pour_recup_motif= element_generation_sequence;
ptr_struct_seq tete_liste_pour_insertion_motif= element_generation_sequence;
ptr_struct_seq element_generation_sequence_next= NULL;
ptr_liste_motif element_motif= malloc(sizeof(liste_chaine_motif));
ptr_liste_motif tete_liste_motif= element_motif;



/* GENERATION ALEATOIRE DE DIX SEQUENCES POUR TEST */
// D�but du code:
printf("aaaaaaargh");
for(i=1; i<=nb_sequence; i++){ //Boucle qui permet de cr�� les 10 s�quences du fichier fasta
	element_generation_sequence->numero_sequence= i;
	for (j=0; j<30; j++){ //boucle qui permet de cr�er les s�quences de 30 nucl�otides de long
            random_nucleotide(&nucleo);
            element_generation_sequence->sequence[j]= nucleo;
		// on utilise la fonction random_nucl�otide pour choisir al�atoirement les 30 nucl�otides.
	}
	element_generation_sequence->sequence[30]= '\0';
	if(i!=10){ // Tant qu'on est pas au dernier �l�ment on continue de cr�er des briques
		element_generation_sequence_next=malloc(sizeof(structure_sequence)); //on cr�e un �l�ment de type structure_s�quence
		element_generation_sequence-> next_sequence= element_generation_sequence_next; // le pointeur de la chaque brique pointe sur l'�l�ment suivant
		element_generation_sequence=element_generation_sequence_next; //on passe de la brique fra�chement cr��e � la suivante
	}else{ //cas du dernier �l�ment de la liste
		element_generation_sequence->next_sequence= NULL;



	}

}
insert_motif(&tete_liste_pour_insertion_motif, nb_sequence, &tete_liste_motif);

generation_masque(longueur_masque, &masque, nb_fenetre);
parcours_masque(longueur_masque, &masque, nb_fenetre, nb_sequence, &tete_liste_pour_parcours_masque, &tete_liste_kmer, &tete_liste_sequence, &tete_liste_pos);
kmer_present_dans_chaque_sequence(nb_sequence, &tete_liste_kmer2, &tete_liste_sequence2, &tete_liste_pos2, &tete_liste_pour_recup_motif, &tete_liste_kmer_selectionne, &tete_liste_motif_PSSM);
affichage_dictionnaire_kmer(&tete_liste_kmer3, &tete_liste_sequence3, &tete_liste_pos3);
affichage_motif_selectionne(&tete_liste_kmer_selectionne2, &tete_liste_motif_PSSM2);
// On calculera la PSSM seulement pour les kmers qui sont présent dans plus de 7 sequences:
while (tete_liste_kmer4 != NULL)
{
	if (tete_liste_kmer4->nb_sequence>=7)
	{
		calcul_PSSM(&tete_liste_kmer_selectionne_pour_calcul, &tete_liste_motif_PSSM_pour_calcul, &file_info, &ptr_matrice_PSSM);
	}
	tete_liste_kmer4= tete_liste_kmer4->suiv_kmer;
}
//printf("\npleaase: %s", p_motif_consensus->motif_consensus);
//calcul_ST1()



return 0;
}
