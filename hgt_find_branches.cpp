//================================================================================================
//=  HGT-DETECTION v3.3.1b
//=  Authors : Alix Boc and Vladimir Makarenkov
//=  Date : June 2009
//=
//=  Description : Cette version permet de generer le fichier outputweb. Les transferts sont déja
//=	 connus.
//=
//================================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>

#pragma warning(disable:4996)

#include "structures.h"
#include "utils_tree.cpp"
#include "fonctions.cpp"

#define binaireSpecies 0 
#define binaireGene    1

void traiterSignal(int sig){
	printf("\nMESSAGE : SEGMENTATION FAULT #%d DETECTED",sig);
	printf("\nUse valgrind or gdb to fix the problem");
	printf("\n");
	exit(-1);
}

//========================================================================================================
//============================================ MAIN ======================================================
//========================================================================================================
int main(int nargc,char **argv){

	struct InputTree SpeciesTree;				    //== initial species tree
	struct InputTree SpeciesTreeCurrent;		//== initial species tree
	struct InputTree FirstTree;
	struct InputTree FirstTree2;
	struct InputTree GeneTree;					    //== initial gene tree
	struct InputTree SpeciesTreeRed;			  //== reduced species tree
	struct InputTree GeneTreeRed;				    //== reduced gene tree
	struct ReduceTrace aMap;					      //== mapping structure between species tree and recuded species tree
	struct InputTree geneTreeSave;
	struct HGT * bestHGTRed = NULL;				  //== list of HGT for the reduced tree
	struct HGT * bestHGT = NULL;				    //== list of HGT for the normale tree
	struct HGT * outHGT = NULL;
	struct HGT * bestHGTmulticheck = NULL;
	int nbHGT_boot;
	int first = 1,k,l;
	int cpt_hgt=0,i,j,tmp,nbTree=0;
	int bootstrap = 0;
	int multigene = 0;
	int nbHgtFound = 0;
	struct CRITERIA * multicheckTab=NULL;
	struct CRITERIA aCrit;						       //== struture of all the criteria
	struct DescTree *DTSpecies,					       //== structure of submatrices for the species tree
					*DTGene;					       //== structure of submatrices for the gene tree
	struct Parameters param;
	FILE *in,*out;
	int max_hgt,nbHGT;
	int ktSpecies;
    int trivial = 1;
	int *speciesLeaves = NULL;
	int RFref;
	int imc;	
	char *mot = (char*)malloc(100);
	int nomorehgt=0;
	initInputTree(&geneTreeSave);


	//== read parameters
	printf("\nhgt : reading options");
	if(readParameters(&param,argv,nargc)==-1){
		printf("\nhgt : no options specified, see the README file for more details\n");
		exit(-1);
	}

	rand_bootstrap = param.rand_bootstrap;
	
	signal(SIGSEGV,traiterSignal);
	
	//== open the input file
	if((in=fopen(param.inputfile,"r"))==NULL){
		printf("\nhgt : The file %s does not exist",param.inputfile);
		exit(-1);
	}
	if(strcmp(param.speciesroot,"file") == 0){
		if(!file_exists(param.speciesRootfileLeaves)){
			printf("\nhgt : The file %s does not exist",param.speciesRootfileLeaves);
			exit(-1);
		}
	}
	if(strcmp(param.generoot,"file") == 0){
		if(!file_exists(param.geneRootfileLeaves)){
			printf("\nhgt : The file %s does not exist",param.geneRootfileLeaves);
			exit(-1);
		}
	}
	if((in=fopen(param.inputfile,"r"))==NULL){
		printf("\nhgt : Cannot open input file (%s)",param.inputfile);
		exit(-1);
	}
	
	//== open the bootstrapFile
	if(strcmp(param.bootstrap,"yes") == 0){
		bootstrap = 1;
	}
	if(strcmp(param.multigene,"yes") == 0){
		multigene = 1;
		if(strcmp(param.speciesroot,"file"))
			strcpy(param.speciesroot,"midpoint");
	}
	remove(param.hgtResultFile);

	initInputTree(&FirstTree);
	initInputTree(&SpeciesTreeCurrent);

	FILE * results, *results_bouba; // = fopen(param.results,"w+");
	
	if((results = fopen(param.results,"w+"))==NULL){
		printf("PROBLEME AVEC RESULTS1");
		exit(0);
	}
	if((results_bouba = fopen(param.results_bouba,"w+"))==NULL){
		printf("PROBLEME AVEC RESULTS2");
		exit(0);
	}

	
//==============================================================================
//============================= LECTURE DES ARBRES =============================
//==============================================================================
	
	printf("\nhgt : reading the input file");
	tmp = readInputFile(in, param.input/*,&SpeciesTree,&GeneTree*/,param.errorFile);

	if(tmp==-1) {
		printf("\nCannot read input data !!\n");
		exit(-1);
	}
  
	
	cpt_hgt = 0;

	initInputTree(&SpeciesTree);
	initInputTree(&GeneTree);
	initInputTree(&SpeciesTreeRed);
	initInputTree(&GeneTreeRed);

	//== lecture des matrices ou chaines newick en entree
	if(readInput(SPECIE,param.input,&SpeciesTree) == -1){ printf("\nError in species tree\n"); exit(-1);}
	if(readInput(GENE,param.input,&GeneTree) == -1){ printf("\nError in gene tree\n"); getchar(); exit(-1);}

	TrierMatrices(GeneTree.Input,GeneTree.SpeciesName,SpeciesTree.SpeciesName,SpeciesTree.size);
	
	NJ(SpeciesTree.Input,SpeciesTree.ADD,SpeciesTree.size);
	NJ(GeneTree.Input,GeneTree.ADD,GeneTree.size);

	//== construction des differentes représentation des arbres (adjacence,aretes,longueur,degre)
	CreateSubStructures(&SpeciesTree,1,binaireSpecies);
	CreateSubStructures(&GeneTree,1,binaireGene);

//for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++)
	//printf("\n%d--%d : %lf",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],SpeciesTree.LONGUEUR[i-1]);

//==============================================================================
//============================ GESTION DES RACINES =============================
//==============================================================================
  //== selection de la racine
	printf("\nhgt : adding the tree roots");
	if(strcmp(param.load,"yes") == 0 ){
		//printf("On charge les fichiers");

		chargerFichier(&SpeciesTree,param.speciesTree,param.speciesRootfile);
		//addRoot(&GeneTree,NULL,GeneBranch,param.generoot,param.geneRootfile,NULL);
		chargerFichier(&GeneTree,param.geneTree,param.geneRootfile);
	}
	else{		
		if((strcmp(param.version,"web")==0) && (strcmp(param.printWeb,"yes")==0)){
		    saveTree(param.speciesTreeWeb,SpeciesTree,bestHGT,0,cpt_hgt,"",param.scenario,NULL);
		    saveTree(param.geneTreeWeb,GeneTree,bestHGT,0,cpt_hgt,"",param.scenario,NULL);
	   }

		if(first == 1){
			//FILE *froot;
			int nbBranche,leave;			
		}
		if(speciesLeaves == NULL){
			speciesLeaves = (int*) malloc(SpeciesTree.size*sizeof(int));
			speciesLeaves[0] = -1;
		}
		
		if(SpeciesTree.Root == -1) addRoot(&SpeciesTree,NULL,SpeciesBranch,param.speciesroot,param.speciesRootfile,param.speciesRootfileLeaves,NULL,param.version);
		if(GeneTree.Root == -1) addRoot(&GeneTree,NULL,GeneBranch,param.generoot,param.geneRootfile,param.geneRootfileLeaves,NULL,param.version);
	
		if(strcmp(param.viewtree,"yes")==0)
			exit(0);
	}

	nbTree++;

	if(SpeciesTree.size > GeneTree.size) max_hgt = 4*GeneTree.size * GeneTree.size;
	else max_hgt = 4*SpeciesTree.size * SpeciesTree.size;

	bestHGTRed = (struct HGT*)malloc(max_hgt*sizeof(struct HGT)); 
	bestHGT = (struct HGT*)malloc(max_hgt*sizeof(struct HGT));

	for(i=0;i<max_hgt;i++){
		bestHGTRed[i].listSource = NULL;
		bestHGTRed[i].listDestination = NULL;
		bestHGT[i].listSource = NULL;
		bestHGT[i].listDestination = NULL;
	}

	if(first==1 && strcmp(param.scenario,"multiple")!=0 && strcmp(param.mode,"multicheck")==0)
		multicheckTab = (struct CRITERIA *)malloc(max_hgt*sizeof(struct CRITERIA)); 

	if(first==1){
		multicheckTab[0].m = 0;
		//imc = 1;
		imc=0;
	}
	InitCriteria(&aCrit,SpeciesTree.size);		

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);

	if(first ==1){
		FirstTree.ADD=NULL;
		FirstTree.ARETE=NULL;
		copyInputTree(&FirstTree,SpeciesTree,1,1);
		AdjustBranchLength(&FirstTree,GeneTree,binaireSpecies,1);
	}	
	
	InitCriteria(&aCrit,SpeciesTree.size);
	computeCriteria(FirstTree.ADD,GeneTree.ADD,FirstTree.size,&aCrit,FirstTree.LONGUEUR,FirstTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	
	AdjustBranchLength(&SpeciesTree,GeneTree,binaireSpecies,1);
	
	if((bootstrap != 1) || (bootstrap==1 && first==1)){
		fprintf(results,"%d,%lf,%lf\n",aCrit.RF,aCrit.LS,aCrit.BD);
		RFref=aCrit.RF;
	}

	/*printf("\n");
	for(i=1;i<=2*SpeciesTree.size-3;i++){
		printf("\n%d-%d --> %lf",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],SpeciesTree.LONGUEUR[i-1]);
	}*/
	
	

	SpeciesTreeCurrent.ADD=NULL;
	SpeciesTreeCurrent.ARETE=NULL;
	copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,0);
	AdjustBranchLength(&SpeciesTreeCurrent,GeneTree,binaireSpecies,1);
	AdjustBranchLength(&SpeciesTree,GeneTree,binaireSpecies,1);

	DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);	
	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);


	printf("\nhgt : pre-treatment process");
	//printf("\nPre-traitement .......");
	//== date : 26 janvier 2008
	//== l'appel a cette fonction permet de recreer le meme sous arbre dans l'arbre de gene et l'arbre d'especes
	//== si les distances entre les especes sont nulles.
	int tousLesCasSontTraitees = TRUE; //FALSE;
	int *tab_tous_les_sommets = (int*)malloc((SpeciesTree.size+1) * sizeof(int));
	int *tab_sommets_selectionnes = (int*)malloc((SpeciesTree.size+1) * sizeof(int));
	int *tab_branches = (int*)malloc( 4*(GeneTree.size+1) * sizeof(int));
	int nb_branches;
	int temoin_nouveau_cas;
	struct HGT aHGT;
	
	for(i=1;i<=SpeciesTree.size;i++) tab_tous_les_sommets[i] = 0;
	tab_tous_les_sommets[0] = FALSE;
	
	while(tousLesCasSontTraitees == FALSE){
		
		ListeSommets_taille_0(GeneTree.Input,tab_tous_les_sommets,GeneTree.size-1);
		tousLesCasSontTraitees = tab_tous_les_sommets[0];
		
		if(tousLesCasSontTraitees == FALSE){
			tab_sommets_selectionnes[0] = 0;
			temoin_nouveau_cas=0;	
			for(i=1;i<GeneTree.size;i++){
				if(tab_tous_les_sommets[i] == 1){
					if(temoin_nouveau_cas == 0){
						temoin_nouveau_cas=1;
						//printf("\nTraitement : ");
					}
					tab_tous_les_sommets[i] = 2;
					tab_sommets_selectionnes[0] = tab_sommets_selectionnes[0] + 1;
					tab_sommets_selectionnes[tab_sommets_selectionnes[0]] = i;
					//printf("%d(%s) ",i,GeneTree.SpeciesName[i]);
				}
			}
			if(tab_sommets_selectionnes[0] > 0){
				ListesBranchesPourHGT(tab_sommets_selectionnes,GeneTree.ARETE,GeneTree.size,DTGene,tab_branches,&nb_branches);
				while(findBestHGT_nombreLimite(DTGene,DTSpecies,tab_branches,nb_branches,GeneTree,SpeciesTree,param,&aHGT) > 0){
					applyHGT(SpeciesTree.ADD,&GeneTree,aHGT.source,aHGT.destination);
					AdjustBranchLength(&GeneTree,SpeciesTree,0,1);
					deleteBipartition(DTGene,GeneTree);
					DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
					RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);	
					
					ListesBranchesPourHGT(tab_sommets_selectionnes,GeneTree.ARETE,GeneTree.size,DTGene,tab_branches,&nb_branches);
				}

			}
			
		}		
	}
	free(tab_tous_les_sommets);
	free(tab_sommets_selectionnes);
	free(tab_branches);
	
	deleteBipartition(DTGene,GeneTree);
	deleteBipartition(DTSpecies,SpeciesTree);
	
	//initInputTree(&GeneTree);
	//copyInputTree(&GeneTree,tmpGeneTree,1,1);
	
	GeneTree.size = GeneTree.size-1;
	CreateSubStructures(&GeneTree,1,binaireGene);
	for(i=1;i<=GeneTree.size+1;i++){
		GeneTree.ADD[GeneTree.size+1][i] = GeneTree.ADD[i][GeneTree.size+1] = epsilon;
	}
	//CreateSubStructures(&GeneTree,1,binaireGene);
	
	AdjustBranchLength(&SpeciesTree,GeneTree,binaireSpecies,1);
	/*
	printf("\n\nSpecies Tree : ");
		for(i=1;i<=SpeciesTree.size;i++){
			printf("\n%s\t",SpeciesTree.SpeciesName[i]);
			for(j=1;j<=SpeciesTree.size;j++){
				printf("%lf ",SpeciesTree.ADD[i][j]);
			}
		}
		printf("\n\nGene Tree : ");
		for(i=1;i<=GeneTree.size;i++){
			printf("\n%s\t",GeneTree.SpeciesName[i]);
			for(j=1;j<=GeneTree.size;j++){
				printf("%lf ",GeneTree.ADD[i][j]);
			}
		}
		
	printf("\n%d",SpeciesTree.size);
	for(i=1;i<=SpeciesTree.size;i++){
		printf("\n%d\t",i);
		for(j=1;j<=SpeciesTree.size;j++){
			printf("%lf ",SpeciesTree.ADD[i][j]);
		}
	}	*/
	DTGene = (struct DescTree*)malloc(3*(2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartitionSansRacine(GeneTree.ARETE,GeneTree.ADD,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);	
		
	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
	
	InitCriteria(&aCrit,SpeciesTree.size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	printf("\nhgt : %d %d %lf\n",aCrit.RF,SpeciesTree.size,aCrit.LS);
	
	
	
	
	
	
	//=====================================================================================================================
	// ./hgt_find_branches -inputfile=input -speciresroot=file(speciesRootLeaves.txt)
	//=====================================================================================================================
	
	for(int i=1; i<=SpeciesTree.size; i++){	
		printf("\n%d-%s",i,SpeciesTree.SpeciesName[i]);
	}
	
	char path[50];
	strcpy (path,param.path);
	FILE * transferts;
	if((transferts = fopen( strcat(path, "outputCons.txt" ),"r") )==NULL){
		printf("PROBLEME AVEC outputCons");
		exit(0);
	}
	strcpy (path,param.path);
	FILE * validation;
	if((validation = fopen( strcat(path, "outputConsValide.txt" ),"w") )==NULL){
		printf("PROBLEME AVEC validation");
		exit(0);
	}

	int val,no,cpt=1;

	multicheckTab[0].m = 1;	

	multicheckTab[1].nbHgtFound = 0;

	char * str_lang = (char*) malloc(100);

	while(fscanf(transferts,"%d",&no) != -1){

		fprintf(validation ,"dans le while et no = %d\n", no);
		if(no > 0){
			cpt_hgt++;
			printf("\nTransfert #%d\nSource (%d-%d) = ",cpt++,no,cpt_hgt);
			bestHGT[cpt_hgt].listSource = (int*)malloc((no+1)*sizeof(int));
			bestHGT[cpt_hgt].listSource[0] = no;
			for(int i=1;i<=no;i++){
				fscanf(transferts,"%s",str_lang);
				printf("%s",str_lang);
				for(int j=1; j<=SpeciesTree.size; j++){	
					if(strcmp(SpeciesTree.SpeciesName[j],str_lang) == 0){
						printf("(%d) ",j);
						bestHGT[cpt_hgt].listSource[i] = j;
					}
				}
			}
			fscanf(transferts,"%d",&no);
			printf("\nDestination(%d-%d) = ",no,cpt_hgt);
			bestHGT[cpt_hgt].listDestination = (int*)malloc((no+1)*sizeof(int));
			bestHGT[cpt_hgt].listDestination[0] = no;
			for(int i=1;i<=no;i++){
				fscanf(transferts,"%s",str_lang);
				printf("%s ",str_lang);
				for(int j=1; j<=SpeciesTree.size; j++){	
					if(strcmp(SpeciesTree.SpeciesName[j],str_lang) == 0){
						printf("(%d) ",j);
						bestHGT[cpt_hgt].listDestination[i] = j;
					}
				}
			}
			fprintf(validation ,"Avant le apply !!!\n");
			val = applyHGTSANSRACINE2(GeneTree.ADD,&SpeciesTree,&bestHGT[cpt_hgt]);
			
			fprintf(validation ,"Apres le apply !!!%d\n",val);
			if (val != -1){
			//fprintf(validation ,"dans le if !!!%d\n",val);
				bestHGT[cpt_hgt].valide = 1;
				multicheckTab[1].nbHgtFound ++;
			}
		}
		else{
			fprintf(validation ,"dans le else !!!%d\n",val);
			fscanf(transferts,"%d",&no);
			printf("\nEtat = %d",no);
		}
	}
	
	fclose(transferts);
	fclose(validation);
	
	
	for(int i=1;i<=cpt_hgt;i++){
		printf("\nTransfert #%d",i);
		printf("\nSource : ");
		for(int j=1; j<=bestHGT[i].listSource[0]; j++){	
			printf("%d ",bestHGT[i].listSource[j]);
		}
		printf("\nDestination : ");
		for(int j=1; j<=bestHGT[i].listDestination[0]; j++){	
			printf("%d ",bestHGT[i].listDestination[j]);
		}
		printf("\n%d--%d -> %d--%d",bestHGT[i].source_A,bestHGT[i].source_B,bestHGT[i].dest_A,bestHGT[i].dest_B);
	}
	

	
	
	printf("\nFirst Tree");
	for(i=1;i<=2*SpeciesTree.size-3;i++){
		printf("\n%d-%d --> %lf",FirstTree.ARETE[2*i-1],FirstTree.ARETE[2*i-2],FirstTree.LONGUEUR[i-1]);
	}
	
//==============================================================================
//============================ TRAITEMENT DES RESULTATS ========================
//==============================================================================
	printf("\nhgt : formatting the results");
	//outHGT = (struct HGT*)malloc(2*param.nbhgt*sizeof(struct HGT)); 
	outHGT = (struct HGT*)malloc(max_hgt*sizeof(struct HGT));
	nbHGT = formatResult(bestHGT,cpt_hgt,outHGT,FirstTree);
	printHGT(results,results_bouba,param.stepbystep,multicheckTab,param.mode,RFref,/*out,*/FirstTree,outHGT,nbHGT,NULL,param.subtree,param.bootmin);

	//if((strcmp(param.version,"web")==0) && (strcmp(param.printWeb,"yes")==0)){
		saveTree(param.outputWeb,FirstTree,outHGT,1,nbHGT,param.subtree,param.scenario,NULL);
	//}
	remove(param.noMoreHgtfile);
	if(nomorehgt == 0){
		FILE * out = fopen(param.noMoreHgtfile,"w+");
		fclose(out);
	}
		
//==============================================================================
//============================ LIBERATION DE LA MEMOIRE ========================
//==============================================================================		
	deleteBipartitionSansRacine(DTGene,GeneTree.size);
	FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
	FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
	FreeMemory_InputTree(&SpeciesTreeCurrent,SpeciesTreeCurrent.size);
	FreeMemory_InputTree(&GeneTree,GeneTree.size);
	FreeMemory_InputTree(&SpeciesTree,SpeciesTree.size);
	if(bootstrap !=1)
		FreeMemory_InputTree(&FirstTree,FirstTree.size);
	FreeCriteria(&aCrit,SpeciesTree.size);
	
	for(i=1;i<=cpt_hgt;i++){
		free(bestHGT[i].listSource);
		free(bestHGT[i].listDestination);
	}
	free(bestHGTRed);
	free(bestHGT); 
	
  fclose(results);
 
 //printf("\nvaleur retournee : %d \n",nbHGT);
  printf("\nhgt : number of HGT(s) found = %d \nhgt : end of computation, check the file results.txt for the program output\n",nbHGT);
  exit(nbHGT);
}
