

//------------------------------------------------------------------------------------------------------------
// lecture du fichier contenant les variables souhaitées, stockage de ces variables pour utilisation
void importer_parametres(int* longueur_masque, int* d, int* nb_fenetre, int* nb_masques)
{
	FILE* ptr_fichier; //creation d'un pointeur sur le fichier
  char nom_fichier[30];

	printf("Quel est le nom du fichier que vous voulez utiliser pour importer les paramètres du programme ?\n");
	scanf("%s", nom_fichier);
  ptr_fichier = fopen(nom_fichier, "r"); //ouverture du fichier

	if( ptr_fichier == NULL) { free(ptr_fichier); return; }// si le fichier est vide on sort de la fonction

  fscanf(ptr_fichier, "#longueur du motif à identifier (IE. taille du masque): %d\n", longueur_masque);
  fscanf(ptr_fichier, "#nombre maximal de substitutions autorisées: %d\n", d);
  fscanf(ptr_fichier, "#nombre de fenêtres dans les masques utilisés: %d\n", nb_fenetre);
  fscanf(ptr_fichier, "#nombre de masques à générer: %d\n", nb_masques);

  fclose(ptr_fichier); //fermeture du fichier
}
//------------------------------------------------------------------------------------------------------------
// Récupère les séquences et les stocke dans une liste chainée
void importer_sequences_fasta( TPtr_info_ensemble_sequences* ptr_info, TPtr_ensemble_sequences* ptr_ensemble )
{
  FILE* ptr_fichier_fasta;
  char nom_fichier_fasta[30];
  char contenu_ligne[128];
  char c;
  int cptr = 0;

	// déclaration des pointeurs permettant de créer la liste chainée de séquences
  TPtr_ensemble_sequences p_new;
  TPtr_ensemble_sequences p = *ptr_ensemble;
  (*ptr_info)->nb_sequences = 0;

	printf("Quel est le nom du fichier que vous voulez utiliser pour importer les séquences ?\n");
  scanf("%s", nom_fichier_fasta);

  ptr_fichier_fasta = fopen(nom_fichier_fasta, "r");
  if( ptr_fichier_fasta == NULL) { free(ptr_fichier_fasta); return; }// si le fichier est vide on sort de la fonction

	do //cette boucle permet de compter les séquences
	{
		c = fgetc (ptr_fichier_fasta); //lecture du fichier caractère par caractère
		if (c == '>') cptr++;  //incrémentation pour chaque nom de séquences
	} while (c != EOF); // l'action est répétée jusqu'à ce que la fin du fichier soit atteinte

	ptr_fichier_fasta = fopen(nom_fichier_fasta, "r");
	do
	{
		fscanf(ptr_fichier_fasta,"%s\n", contenu_ligne); // lecture du fichier ligne par ligne
	  if ( contenu_ligne[0] == '>')
	  {
			if (((*ptr_info)->nb_sequences) > 0)
			{
				p_new = malloc ( sizeof(TEnsemble_Sequences));
				p->suiv_seq = p_new;
				p = p_new;
				if ( (*ptr_info)->nb_sequences  == cptr) { p->suiv_seq = NULL;}
			}
			(*ptr_info)->nb_sequences += 1;
			strcpy(p->nom_seq, contenu_ligne);
		}else{
			if ( p->sequence == NULL){
				strcpy(p->sequence, contenu_ligne);
	    }else{
	      strcat(p->sequence, contenu_ligne);
	    }
	  }
	}	while(!feof(ptr_fichier_fasta));
  fclose(ptr_fichier_fasta);
}
//------------------------------------------------------------------------------------------------------------
void afficher_sequences(TPtr_info_ensemble_sequences* ptr_info, TPtr_ensemble_sequences* ptr_ensemble )
{
	FILE * ptr_fichier;
  TPtr_ensemble_sequences p = *ptr_ensemble;
	// TPtr_ensemble_sequences p_prec = NULL;
	int cptr = 0;
	printf("fonction afficher_sequences\n");
	printf("%d \n", (*ptr_info)->nb_sequences);

	ptr_fichier = fopen ( "verif_dico_fasta.txt" ,"w");
	while ( cptr < ((*ptr_info)->nb_sequences))
	{
		fputs( p->nom_seq, ptr_fichier);
		fputs( "\n" ,ptr_fichier);
		fputs( p->sequence,ptr_fichier);
		fputs( "\n" ,ptr_fichier);
		printf("nom_seq: %s\n", p->nom_seq );
		printf("%s\n", p->sequence);
		p = p->suiv_seq;
		cptr ++;
	}
}
