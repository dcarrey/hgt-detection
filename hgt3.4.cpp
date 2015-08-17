//================================================================================================
//=  HGT-DETECTION v3.4
//=  Authors : Alix Boc and Vladimir Makarenkov
//=  Date : March 2012
//=
//=  Description : This program detect horizontal gene transfer (HGT). As input it takes 2 
//=  trees: a species tree and a gene tree. the goal is to transform the species tree
//=  into the gene tree following a transfer scenario. There are 3 criteria : the robinson and
//=  Foulds distance, the least-square criterion and the bipartition distance. We also use the
//=  subtree constraint. With this version we can now perform simulation.
//=
//=	 input   : file with species tree and gene tree in the newick format.
//=			   In case of simulation, the species tree and all the gene trees in the same file in
//=            the phylip format or newick string
//=  output  : a list of HGT and the criteria values for each one.
//=	 options :
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
	struct InputTree FirstTreeB;
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
	int cpt_hgt,i,j,tmp,nbTree=0;
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
		if(!file_exists(param.speciesRootfileLeaves) && !file_exists(param.speciesRootfile)){
			printf("\nhgt : The file %s does not exist",param.speciesRootfileLeaves);
			exit(-1);
		}
	}
	if(strcmp(param.generoot,"file") == 0){
		if(!file_exists(param.geneRootfileLeaves) && !file_exists(param.geneRootfile)){
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
		printf("PROBLEME AVEC RESULTS");
		exit(0);
	}
	if((results_bouba = fopen(param.results_bouba,"w+"))==NULL){
		printf("PROBLEME AVEC RESULTS");
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

	InitCriteria(&aCrit,SpeciesTree.size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	printf("%d,%lf,%lf\n",aCrit.RF,aCrit.LS,aCrit.BD);

	//== construction des differentes représentation des arbres (adjacence,aretes,longueur,degre)
	CreateSubStructures(&SpeciesTree,1,binaireSpecies);
	CreateSubStructures(&GeneTree,1,binaireGene);

	InitCriteria(&aCrit,SpeciesTree.size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	printf("%d,%lf,%lf\n",aCrit.RF,aCrit.LS,aCrit.BD);
	
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
	else
	{		
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
		bestHGTRed[i].valide = 1;
		bestHGT[i].listSource = NULL;
		bestHGT[i].listDestination = NULL;
	}

	/*if(first==1 && strcmp(param.scenario,"multiple")!=0 && strcmp(param.mode,"multicheck")==0)
		multicheckTab = (struct CRITERIA *)malloc(max_hgt*sizeof(struct CRITERIA)); 

	if(first==1){
		multicheckTab[0].m = 0;
		//imc = 1;
		imc=0;
	}
	*/
	
	if(first==1)
	{
		multicheckTab = (struct CRITERIA *)malloc(max_hgt*sizeof(struct CRITERIA)); 
	
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

	
	//====================================================================
	//== Ajout de transferts avant processus de detection
	//====================================================================
	//printf("\nhgt : stepbystep=%s",param.stepbystep);
	//if(strcmp(param.stepbystep,"yes") == 0){
	//== si le fichier prehgtfile existe, 
	if(file_exists(param.prehgtfile))
	{
		FILE *prehgt = fopen(param.prehgtfile,"r");
		
		printf("\nhgt : adding pre-hgt");
		int pos_source,pos_dest,step,valide,newStep=-1;
		int a1,a2,b1,b2;
		int dedans=0,dedans2=0;
		int nbHGTadd=0;
		while(fscanf(prehgt,"%d%d%d%d%d%d",&step,&a1,&a2,&b1,&b2,&valide) != -1){
			dedans=0;
			SpeciesTreeCurrent.ADD=NULL;
			SpeciesTreeCurrent.ARETE=NULL;
			copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,1);
			
			DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
			RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);

			
			//printf("\nhgt : %d %d %d %d %d %d",step,a1,a2,b1,b2,valide);
			
			if (valide == 2){
				int tmp=a1; a1=b1; b1=tmp;
				tmp=a2; a2=b2; b2=tmp;
				//printf("\nhgt : %d %d %d %d %d %d",step,a1,a2,b1,b2,valide);
				valide=6;
			}
			pos_source = pos_dest = -1;
			//AdjustBranchLength(&SpeciesTree,GeneTree,binaireSpecies,1);
			for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++){
					//printf("\n%d--%d : %lf",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],SpeciesTree.LONGUEUR[i-1]);
					if((SpeciesTree.ARETE[2*i-1] == a1 && SpeciesTree.ARETE[2*i-2] == a2) || (SpeciesTree.ARETE[2*i-2] == a1 && SpeciesTree.ARETE[2*i-1] == a2))
						pos_source = i;
					if((SpeciesTree.ARETE[2*i-1] == b1 && SpeciesTree.ARETE[2*i-2] == b2) || (SpeciesTree.ARETE[2*i-2] == b1 && SpeciesTree.ARETE[2*i-1] == b2))
						pos_dest = i;
			}
			printf("\nhgt : (branches) %d %d",pos_source,pos_dest);
			if(pos_source != -1 && pos_dest != -1){
				if(isAValidHGT(SpeciesTree,pos_source,pos_dest) == 1){
					dedans=1;
					if(step != newStep){
						if(newStep >= 0){
							imc++;
						}
						else{
							imc=1;
						}
						newStep = step;
						multicheckTab[0].m ++;	
						multicheckTab[imc].nbHgtFound =0;
					}
					
					//nbHGTadd++;
					printf("\nhgt : (branches) %d %d",pos_source,pos_dest);
					applyHGT(NULL,&SpeciesTree,pos_source,pos_dest);
					AdjustBranchLength(&SpeciesTree,GeneTree,binaireSpecies,1);
					computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
					printf("\nhgt : %d-%d -> %d-%d",a1,a2,b1,b2);
					printf("\nhgt : CRITERIA RF=%d,LS=%lf,BD=%lf\n",aCrit.RF,aCrit.LS,aCrit.BD);
					
					multicheckTab[imc].nbHgtFound ++;
					cpt_hgt ++;
					bestHGT[cpt_hgt].source_A = a1;
					bestHGT[cpt_hgt].source_B = a2;
					bestHGT[cpt_hgt].dest_A = b1;
					bestHGT[cpt_hgt].dest_B = b2;
					bestHGT[cpt_hgt].valide = valide;
					bestHGT[cpt_hgt].crit = aCrit;
					bestHGT[cpt_hgt].source = pos_source;
					bestHGT[cpt_hgt].destination = pos_dest;
					bestHGT[cpt_hgt].sequence = imc;
					if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A) || 
						(bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B) || 
						(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A) || 
						(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B)){
						bestHGT[cpt_hgt].trivial = 1;
						bestHGT[cpt_hgt].valide = TRIVIAL;
					}
					if(valide == 3){
						bestHGT[cpt_hgt].valide = valide;
					}
					
					multicheckTab[imc].LS = aCrit.LS;
					multicheckTab[imc].RF = aCrit.RF;
					multicheckTab[imc].BD = aCrit.BD;
					multicheckTab[imc].rLS = -1;
					multicheckTab[imc].rRF = -1;
					multicheckTab[imc].rBD = -1;
					findListSpecies(&bestHGT[cpt_hgt],DTSpecies,SpeciesTreeCurrent);
					
					InitCriteria(&aCrit,SpeciesTree.size);
				}
			}
		}
		//if(dedans == 1)
			//imc++;
	}

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
	int tousLesCasSontTraitees = FALSE;
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
	
	InitCriteria(&aCrit,SpeciesTree.size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	//printf("\nhgt : %d %d %lf\n",aCrit.RF,SpeciesTree.size,(double)aCrit.RF/(2*(double)SpeciesTree.size-6));
	
	ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);

	tmp = 0;


	//=============================== DETECTION DES TRANSFERTS ====================================
	if(strcmp(param.version,"consol")==0)
	{
		//printf("\nDetection in progress",nbTree);		
		printf("\n===============================================");
		printf("\n| CRITERIA VALUES BEFORE DETECTION ");	
		printf("\n| RF distance  = %2d",aCrit.RF);
		printf("\n| LS criterion = %2.1lf",aCrit.LS);
		printf("\n| BD criterion = %2.1lf",aCrit.BD);
		printf("\n===============================================\n");
	}

	if(strcmp(param.scenario,"multiple")==0)
	{
		cpt_hgt = findAllHGT(SpeciesTreeRed,GeneTreeRed,param,bestHGTRed);
		for(i=1;i<=cpt_hgt;i++)
		{
			expandBestHGT(bestHGTRed[i],&bestHGT[i],aMap,DTSpecies,SpeciesTree);
		}
		sortHGT(bestHGT,cpt_hgt,param);
		printf("\nhgt : scenario multiple");
		if(cpt_hgt > param.nbhgt) cpt_hgt = param.nbhgt;
	}
	else if(strcmp(param.scenario,"unique")==0)
	{
		printf("\nhgt : scenario unique");
		if( strcmp(param.subtree,"yes") == 0 && strcmp(param.mode,"multicheck")==0 )
		{
      
			//==============================================================================
			//============ RECHERCHE DE TRANSFERTS - PLUSIEURS PAR TOUR ====================
			//==============================================================================
			printf("\nhgt : start of detection");
			trivial = (SpeciesTree.kt == 0)?0:1;	
			while( findBestHGTtab(SpeciesTreeRed,GeneTreeRed,param,bestHGTRed,&nbHgtFound,&trivial,bootstrap) > 0)
			{
				imc++;
				trivial = (SpeciesTree.kt == 0)?0:1;

				printf("\n\n[%d HGT%s]", nbHgtFound,(nbHgtFound > 1)?"s":"");	
				
				if(first==1)
				{
					multicheckTab[0].m ++; //= nombre d'occurences
					multicheckTab[imc].nbHgtFound = nbHgtFound;
				}
				
				int temoin_zero=-1;
				int cpt_hgt2 = cpt_hgt;
				
				trier(bestHGTRed,nbHgtFound);
				
				for(i=0;i<nbHgtFound;i++)
				{
					if(bestHGTRed[i].crit.RF == 0)
					{
						temoin_zero = i;
					}
				}
				printf("\nHGT-DETECTION : cpt_hgt= %d",cpt_hgt);
				initInputTree(&FirstTree2);
				copyInputTree(&FirstTree2,SpeciesTree,0,0);
				for(i=0;i<nbHgtFound;i++)
				{
					cpt_hgt++;		
					
					expandBestHGT(bestHGTRed[i],&bestHGT[cpt_hgt],aMap,DTSpecies,SpeciesTree);
					printf("\nhgt_red : HGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGTRed[i].source_A,bestHGTRed[i].source_B,bestHGTRed[i].dest_A,bestHGTRed[i].dest_B);
						
					bestHGT[cpt_hgt].sequence = imc;
					bestHGT[cpt_hgt].trivial = 0;
					bestHGTRed[i].listSource = NULL;
					bestHGTRed[i].listDestination = NULL;
					if((temoin_zero != -1)&&(i!=temoin_zero))
					{
						bestHGT[cpt_hgt].valide = 0; 
						printf("\nHGT-DETECTION : un des transfert met RF=0");
						continue;  
					}
					  
					if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A) || 
						(bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B) || 
						(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A) || 
						(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B))
					{
						bestHGT[cpt_hgt].trivial = 1;
						bestHGT[cpt_hgt].valide = TRIVIAL;
					}
									
					if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) )
					{
					   

						bestHGT[cpt_hgt].valide = 0;
						printf("\nhgt : useless=> HGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);

						continue;
					}	
					applyHGT2(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
					// if(bestHGT[cpt_hgt].valide > 0){
						// printf("\nhgt : HGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
						// applyHGT2(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
					// }
				}
				
				// initInputTree(&FirstTree2);
				// copyInputTree(&FirstTree2,SpeciesTree,0,0);
				cpt_hgt -= nbHgtFound;
				printf("\nHGT-DETECTION : cpt_hgt= %d",cpt_hgt);
				for(i=0;i<nbHgtFound;i++)
				{
					cpt_hgt++;
					if(bestHGT[cpt_hgt].valide > 0)
					{
						SpeciesTree.ADD = NULL;
						SpeciesTree.ARETE = NULL;
						copyInputTree(&SpeciesTree,FirstTree2,0,0);
						applyHGT2(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
						
						computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
						bestHGT[cpt_hgt].crit.LS = aCrit.LS;
						bestHGT[cpt_hgt].crit.RF = aCrit.RF;
						bestHGT[cpt_hgt].crit.BD = aCrit.BD;
						printf("\nRF = %d | LS = %1.2lf | BD = %1.2lf",bestHGT[cpt_hgt].crit.RF,bestHGT[cpt_hgt].crit.LS,bestHGT[cpt_hgt].crit.BD);	
						
						SpeciesTree.ADD = NULL;
						SpeciesTree.ARETE = NULL;
						copyInputTree(&SpeciesTree,FirstTree2,0,0);
						applyHGT2(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].destination,bestHGT[cpt_hgt].source);
					
						computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
						bestHGT[cpt_hgt].crit.rLS = aCrit.LS;
						bestHGT[cpt_hgt].crit.rRF = aCrit.RF;
						bestHGT[cpt_hgt].crit.rBD = aCrit.BD;
						printf("\nrRF = %d | rLS = %1.2lf | rBD = %1.2lf",bestHGT[cpt_hgt].crit.rRF,bestHGT[cpt_hgt].crit.rLS,bestHGT[cpt_hgt].crit.rBD);	
					}
				}
				
				SpeciesTree.ADD = NULL;
				SpeciesTree.ARETE = NULL;
				copyInputTree(&SpeciesTree,FirstTree2,0,0);
				cpt_hgt -= nbHgtFound;
				printf("\nHGT-DETECTION : cpt_hgt= %d",cpt_hgt);
				for(i=0;i<nbHgtFound;i++)
				{
					cpt_hgt++;
					if(bestHGT[cpt_hgt].valide != 0){
						printf("\nhgt : HGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
						applyHGT2(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
					}
				}
				
				computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
						
				printf("\n\nCriteria values after this step :");
				printf("\nRF = %d | LS = %1.2lf | BD = %1.2lf\n",aCrit.RF,aCrit.LS,aCrit.BD);	 

				if(first==1)
				{
					multicheckTab[imc].LS = aCrit.LS;
					multicheckTab[imc].RF = aCrit.RF;
					multicheckTab[imc].BD = aCrit.BD;
					multicheckTab[imc].QD = aCrit.QD;
					//imc++;
				}	
				
				if((cpt_hgt >= param.nbhgt)||(aCrit.RF == 0)) break;
						 
				deleteBipartition(DTSpecies,SpeciesTreeCurrent);
				copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,1);
				
				DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
				RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
				FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
				FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
				initInputTree(&SpeciesTreeRed);
				initInputTree(&GeneTreeRed);
				free(aMap.map);
				free(aMap.gene);
				free(aMap.species);
				ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);
				
				if(strcmp(param.stepbystep,"yes") == 0) {
					nomorehgt=1;
					break;
				}	
			}
			free(aMap.map);
			free(aMap.gene);
			free(aMap.species);
			deleteBipartition(DTSpecies,SpeciesTreeCurrent);

			int retour=0;
			int cpt,nbTours,j;
			//printf("\nhgt : multicheckTab[3].nbHgtFound = %d",multicheckTab[3].nbHgtFound);
			printf("\nhgt : Looking for idling hgt");
			retour = DeleteUseLessHGT(cpt_hgt,bestHGT,SpeciesTree,FirstTree); 
			printf("(%d)",retour);
			cpt=1;
			if(first==1)
			{
				printf("\nhgt : multicheckTab[0].m=%d : ",multicheckTab[0].m);
				for(i=1;i<=multicheckTab[0].m;i++){
					nbTours = multicheckTab[i].nbHgtFound;
					printf("(%d) ",nbTours);
					for(j=1;j<=nbTours;j++){
						printf("%d ",cpt);
						if(bestHGT[cpt].valide == 0){
							multicheckTab[i].nbHgtFound --;
						}
						cpt++;
					}
				}
			}
			//printf("\nhgt : recherche des transferts inutiles");
		
		}	
		else
		{
			int initial;
			if(strcmp(param.stepbystep,"no") == 0) 
			{
				initial=1;
					
				while(findBestHGT(initial,SpeciesTreeRed,GeneTreeRed,param,&bestHGTRed[cpt_hgt+1]) > 0)
				{
					initInputTree(&FirstTree2);
					copyInputTree(&FirstTree2,SpeciesTree,0,0);
					cpt_hgt++;	
					
					imc++;
					multicheckTab[0].m ++; //= nombre d'occurences
					multicheckTab[imc].nbHgtFound = 1;
					
					expandBestHGT(bestHGTRed[cpt_hgt],&bestHGT[cpt_hgt],aMap,DTSpecies,SpeciesTree);
					bestHGT[cpt_hgt].trivial = 0;
					if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A) || 
						(bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B) || 
						(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A) || 
						(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B))	
					{
						printf("\n je rentre 1");
						bestHGT[cpt_hgt].trivial = 1;
						bestHGT[cpt_hgt].valide = TRIVIAL;
					}	
					
					printf("\n valeur trivail dans le while1 %d, i = %d",bestHGT[cpt_hgt].valide, cpt_hgt);
					bestHGTRed[cpt_hgt].listSource = NULL;
					bestHGTRed[cpt_hgt].listDestination = NULL;
					
					applyHGT(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
					
					
					initInputTree(&FirstTreeB);
					copyInputTree(&FirstTreeB,SpeciesTree,0,0);
					if(bestHGT[cpt_hgt].valide > 0)
					{
						AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
						computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
						bestHGT[cpt_hgt].crit.LS = aCrit.LS;
						bestHGT[cpt_hgt].crit.RF = aCrit.RF;
						bestHGT[cpt_hgt].crit.BD = aCrit.BD;
						printf("\nRF = %d | LS = %1.2lf | BD = %1.2lf",bestHGT[cpt_hgt].crit.RF,bestHGT[cpt_hgt].crit.LS,bestHGT[cpt_hgt].crit.BD);	
						
						SpeciesTree.ADD = NULL;
						SpeciesTree.ARETE = NULL;
						copyInputTree(&SpeciesTree,FirstTree2,0,0);
						applyHGT(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].destination,bestHGT[cpt_hgt].source);
						AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
						computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
						bestHGT[cpt_hgt].crit.rLS = aCrit.LS;
						bestHGT[cpt_hgt].crit.rRF = aCrit.RF;
						bestHGT[cpt_hgt].crit.rBD = aCrit.BD;
						printf("\nrRF = %d | rLS = %1.2lf | rBD = %1.2lf",bestHGT[cpt_hgt].crit.rRF,bestHGT[cpt_hgt].crit.rLS,bestHGT[cpt_hgt].crit.rBD);	
					}
					
					copyInputTree(&SpeciesTree,FirstTreeB,0,0);
					AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
					computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
					int val = bestHGT[cpt_hgt].valide;
					
					loadCriteria(aCrit,&(bestHGT[cpt_hgt]));
					bestHGT[cpt_hgt].valide = val;
					printf("\n\nCriteria values after this step :");
					printf("\nRF = %d | LS = %1.2lf | BD = %1.2lf\n",aCrit.RF,aCrit.LS,aCrit.BD);	 

					multicheckTab[imc].LS = aCrit.LS;
					multicheckTab[imc].RF = aCrit.RF;
					multicheckTab[imc].BD = aCrit.BD;
					multicheckTab[imc].QD = aCrit.QD;
					
					if(strcmp(param.version,"consol")==0)
					{
						printf("\nHGT #%d %d--%d -> %d--%d",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
						printf("\nRF = %d, LS = %lf, BD = %lf \n",aCrit.RF,aCrit.LS,aCrit.BD);	 
					} 
					if(bestHGT[cpt_hgt].crit.RF == 0) break;// || bestHGT[cpt_hgt].crit.LS < epsilon) break;

					if(cpt_hgt >= param.nbhgt) break;
				
					deleteBipartition(DTSpecies,SpeciesTreeCurrent);
					copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,1);
					DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
					RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
					
					free(aMap.map);
					free(aMap.gene);
					free(aMap.species);
					FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
					FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
					initInputTree(&SpeciesTreeRed);
					initInputTree(&GeneTreeRed);
						
					ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);

				}
			}
			
			initial=0; 
			while( findBestHGT(initial,SpeciesTreeRed,GeneTreeRed,param,&bestHGTRed[cpt_hgt+1]) > 0)
			{	
				initInputTree(&FirstTree2);
				copyInputTree(&FirstTree2,SpeciesTree,0,0);
				cpt_hgt++;
				imc++;
				multicheckTab[0].m ++; //= nombre d'occurences
				multicheckTab[imc].nbHgtFound = 1;
				
				expandBestHGT(bestHGTRed[cpt_hgt],&bestHGT[cpt_hgt],aMap,DTSpecies,SpeciesTree);
              	bestHGT[cpt_hgt].trivial = 0;
          		if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A) || 
             		(bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B) || 
             		(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A) || 
             		(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B))
				{
					printf("\n je rentre 2");
              		bestHGT[cpt_hgt].trivial = 1;
					bestHGT[cpt_hgt].valide = TRIVIAL;
          		}	
				
				printf("\n valeur trivail dans le while2 %d, i = %d",bestHGT[cpt_hgt].valide, cpt_hgt);
				bestHGTRed[cpt_hgt].listSource = NULL;
				bestHGTRed[cpt_hgt].listDestination = NULL;

				applyHGT(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
				initInputTree(&FirstTreeB);
				copyInputTree(&FirstTreeB,SpeciesTree,0,0);
				if(bestHGT[cpt_hgt].valide > 0)
				{
					AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
					computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
					bestHGT[cpt_hgt].crit.LS = aCrit.LS;
					bestHGT[cpt_hgt].crit.RF = aCrit.RF;
					bestHGT[cpt_hgt].crit.BD = aCrit.BD;
					printf("\nRF = %d | LS = %1.2lf | BD = %1.2lf",bestHGT[cpt_hgt].crit.RF,bestHGT[cpt_hgt].crit.LS,bestHGT[cpt_hgt].crit.BD);	
					
					SpeciesTree.ADD = NULL;
					SpeciesTree.ARETE = NULL;
					copyInputTree(&SpeciesTree,FirstTree2,0,0);
					applyHGT(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].destination,bestHGT[cpt_hgt].source);
					AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
					computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
					bestHGT[cpt_hgt].crit.rLS = aCrit.LS;
					bestHGT[cpt_hgt].crit.rRF = aCrit.RF;
					bestHGT[cpt_hgt].crit.rBD = aCrit.BD;
					printf("\nrRF = %d | rLS = %1.2lf | rBD = %1.2lf",bestHGT[cpt_hgt].crit.rRF,bestHGT[cpt_hgt].crit.rLS,bestHGT[cpt_hgt].crit.rBD);	
				}
				
				copyInputTree(&SpeciesTree,FirstTreeB,0,0);
				AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
				computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
				int val = bestHGT[cpt_hgt].valide;
				loadCriteria(aCrit,&(bestHGT[cpt_hgt]));
				bestHGT[cpt_hgt].valide = val;
				printf("\n\nCriteria values after this step :");
				printf("\nRF = %d | LS = %1.2lf | BD = %1.2lf\n",aCrit.RF,aCrit.LS,aCrit.BD);	 

				multicheckTab[imc].LS = aCrit.LS;
				multicheckTab[imc].RF = aCrit.RF;
				multicheckTab[imc].BD = aCrit.BD;
				multicheckTab[imc].QD = aCrit.QD;
					
				
				if(strcmp(param.version,"consol")==0){
					printf("\nHGT #%d %d--%d -> %d--%d",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
					printf("\nRF = %d, LS = %lf, BD = %lf\n",aCrit.RF,aCrit.LS,aCrit.BD);	 
				}
				
				if(cpt_hgt >= param.nbhgt) break;
			
				deleteBipartition(DTSpecies,SpeciesTreeCurrent);
				copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,1);
				free(aMap.map);
				free(aMap.gene);
				free(aMap.species);
				DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
				RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
				FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
				FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
				initInputTree(&SpeciesTreeRed);
				initInputTree(&GeneTreeRed);
			
				ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);
				//printf("==>%d",SpeciesTreeRed.size);
				//printf("<br> prochlorococcus = %d", aMap.map[6]);
				if(strcmp(param.stepbystep,"yes") == 0) {
					nomorehgt=1;
					break;
				}	
			}
			
			first=0;
			printf("\nhgt_valide_mono1="); 
			for(i=1;i<=cpt_hgt;i++)
			{
				if((bestHGT[i].valide >=1)/* && (bestHGT[i].trivial==0)*/)
				{
					if(first==1) printf(","); 
					first=1;	
					printf("%d",bestHGT[i].valide);
				}
			}
			
			free(aMap.map);
			free(aMap.gene);
			free(aMap.species);
			deleteBipartition(DTSpecies,SpeciesTreeCurrent);
			//FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
			//FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
			int retour=0;
			int cpt,nbTours,j;
			//printf("\nhgt : multicheckTab[3].nbHgtFound = %d",multicheckTab[3].nbHgtFound);
			printf("\nhgt : Looking for idling hgt");
			retour = DeleteUseLessHGT(cpt_hgt,bestHGT,SpeciesTree,FirstTree); 
			printf("(%d)",retour);
			cpt=1;
			if(first==1)
			{
				printf("\nhgt : multicheckTab[0].m=%d : ",multicheckTab[0].m);
				for(i=1;i<=multicheckTab[0].m;i++)
				{
					nbTours = multicheckTab[i].nbHgtFound;
					printf("(%d) ",nbTours);
					for(j=1;j<=nbTours;j++)
					{
						printf("%d ",cpt);
						if(bestHGT[cpt].valide == 0)
						{
							multicheckTab[i].nbHgtFound --;
						}
						cpt++;
					}
				}
			}
			
			first=0;
			printf("\nhgt_valide_mono="); 
			for(i=1;i<=cpt_hgt;i++)
			{
				if((bestHGT[i].valide >=1)/* && (bestHGT[i].trivial==0)*/)
				{
					if(first==1) printf(","); 
					first=1;	
					printf("%d",bestHGT[i].valide);
				}
			}
		}
		
		/*printf("\nnbTour=%d",multicheckTab[0].m);
		for(i=1;i<=multicheckTab[0].m;i++){
			printf("\ni=%d,nbHgtFound=%d",i,multicheckTab[i].nbHgtFound);
		}*/
		
				
	}
	
	//==============================================================================
	//============================ TRAITEMENT DES RESULTATS ========================
	//==============================================================================
	printf("\nhgt : formatting the results");
	outHGT = (struct HGT*)malloc(2*param.nbhgt*sizeof(struct HGT)); 
	nbHGT = formatResult(bestHGT,cpt_hgt,outHGT,FirstTree);
	printHGT(results,results_bouba,param.stepbystep,multicheckTab,param.mode,RFref,/*out,*/FirstTree,outHGT,nbHGT,NULL,param.subtree,param.bootmin);
	
	first=0;
	printf("\nhgt_valide="); 
	for(i=1;i<=cpt_hgt;i++)
	{
		if((outHGT[i].valide >=1)/* && (bestHGT[i].trivial==0)*/)
		{
			if(first==1) printf(","); 
			first=1;	
			printf("%d",outHGT[i].valide);
		}
	}
	
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
	deleteBipartition(DTGene,GeneTree);
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
