#!/usr/bin/perl

use strict;
use warnings;
use File::Copy;


#======================================================
#= VERIFICATION DES ARGUMENTS DE LA LIGNE DE COMMANDE
#======================================================
if( scalar @ARGV < 0){
    print "\nErreur\nusage : $0";
    exit 0;
}

my $cmd = "hgt3.4 ";
my $inputfile = "";
my $outputfile="output.txt";
my $bootstrap = "no";
my $consensus = "no";
my $cons=50;
my $path = "./";
my $viewtree="no";
my @tmp_tab;
my @tmp_tab_init;
my %hgt;
my $ligne;
my $nbLines = 5;
my %hgt_number_tab;
my %hgt_description_tab;
my %hgt_compteur_tab;
my %hgt_compteur_tab_reverse;
my %hgt_criterion_tab;
my %hgt_criterion_tab2;
my %hgt_nbHGT_tab;
my @hgt_pos;
my @hgt_pos2;
my $mode;
my $total_hgt;
my $total_trivial;
my $val_retour=0;   #= nombre de hgt trouve
my @hgt_tab;
my $rand_bootstrap = 0;
my $speciesroot = "midpoint";
my $generoot = "midpoint";
my $stepbystep = "no";

#==== READ PARAMETERS ====
foreach my $elt (@ARGV){
  $cmd .= $elt . " " if(($elt !~ "consensus")&&($elt !~ "consvalue"));
  if($elt =~ "bootstrap"){
    @tmp_tab = split("=",$elt);
    $bootstrap = $tmp_tab[1];
    chomp($bootstrap);
    #print STDOUT "$bootstrap";
  }
  
  
   if($elt =~ "consensus"){
    @tmp_tab = split("=",$elt);
    $consensus = $tmp_tab[1];
    chomp($consensus);

  }
  
  if($elt =~ "consvalue"){
    @tmp_tab = split("=",$elt);
    $cons = $tmp_tab[1];
    chomp($cons);
    #print STDOUT "$bootstrap";
  }
  
  if($elt =~ "speciesroot"){
    @tmp_tab = split("=",$elt);
    $speciesroot = $tmp_tab[1];
    chomp($speciesroot);
  }
  if($elt =~ "generoot"){
    @tmp_tab = split("=",$elt);
    $generoot = $tmp_tab[1];
    chomp($generoot);
  }
  if($elt =~ "inputfile"){
    @tmp_tab = split("=",$elt);
    $inputfile = $tmp_tab[1];
  }
  if($elt =~ "path"){
    @tmp_tab = split("=",$elt);
    $path = $tmp_tab[1];
  }
  if($elt =~ "viewtree"){
    @tmp_tab = split("=",$elt);
    $viewtree = $tmp_tab[1];
  }
  if($elt =~ "outputfile"){
    @tmp_tab = split("=",$elt);
    $outputfile = $tmp_tab[1];
  }
  if($elt =~ "help"){
	print_description();
	print_help();
	exit;
  }
  if($elt =~ "stepbystep"){
	@tmp_tab = split("=",$elt);
    $stepbystep = $tmp_tab[1];
  }
}

my $inputfile_bkp = $inputfile;
$inputfile    = "$path" . "$inputfile";
$outputfile = "$path" . "$outputfile";
$cmd          = "./" . $cmd;
#$cmd          = "usagers/" . $cmd;
my $results   = "$path" . "results.txt";
my $firstResults   = "$path" . "firstResults.txt";
my $results2   = "$path" . "results2.txt";
my $firstResults2   = "$path" . "firstResults2.txt";
my $tmp_input = "$path" . "tmp_input.txt";
my $input_no_space    = "$path" . "input_no_space.txt";
my $return_file = "$path" . "return.txt";
my $log_file = "$path" . "log.txt";
my $log_file1 = "$path" . "log1.txt";
my $output_tmp;
my $outputWeb = "$path" . "outputWeb.txt";
my $firstOutputWeb = "$path" . "firstOutputWeb.txt";
my $outputConsensus = "$path" . "outputCons.txt";
my $outputWebConsensus = "$path" . "outputWebCons.txt";
my $generootfile = "$path" . "geneRootLeaves.txt";
my $generootfileBranches = "$path" . "geneRoot.txt";
my $speciesrootfile = "$path" . "speciesRootLeaves.txt";
my $speciesrootfileBranches = "$path" . "speciesRoot.txt";
my $generootfiletmp = "$path" . "geneRootTmp.txt";
my $speciesrootfiletmp = "$path" . "speciesRootTmp.txt";
my $inputfileformated = "$path" . "inputfileformated.txt";
my $prehgtfile = "$path" . "prehgt.txt";
my $hgt = "";
my $premier = 1;
my $nb = 1;
my	$nbHGT=0;
my 	$nbHGT2=0;

 #===== PRINT HEADER =====
  print_title();

  
#= linux like  
`rm -rf $results $outputfile $log_file $return_file`;

#= windows
#`del $results $outputfile $log_file $return_file $outputWeb`;

#===== CHECKING FILES =====
if( $inputfile eq ""){
	print STDOUT "\n\nRUN_HGT : There is no input file";
	exit -1;
}
if( ! -e $inputfile){
	print STDOUT "\n\nRUN_HGT : $inputfile doesn't exist";
	exit -1;
}

if( ($speciesroot eq "file") && ( ! -e $speciesrootfile && ! -e $speciesrootfileBranches) ){
	print STDOUT "\n\nRUN_HGT : $speciesrootfile doesn't exist";
	exit -1;
}

if( ($generoot eq "file") && ( ! -e $generootfile && ! -e $generootfileBranches) ){
	print STDOUT "\n\nRUN_HGT : $generootfile doesn't exist";
	exit -1;
}
   
  #=== LECTURE DE L'ARBRE D'ESPECES ===
  open(IN,"$inputfile") ||  die "Cannot open $inputfile";
  open(OUT,">$inputfileformated") ||  die "Cannot open $inputfileformated";
  while($ligne = <IN>){
	chomp($ligne);
	$ligne =~ s/;/;\n/g;
	if($ligne ne ""){
		print OUT $ligne;
	}
  }
  close(IN);
  close(OUT); 
  open(IN,"$inputfileformated") ||  die "Cannot open $inputfileformated";
  my @trees_tab = <IN>;
  close(IN);
  
  #===========================================================================
  #======================== EXECUTION DU PROGRAMME ===========================
  #===========================================================================
  $cmd .= "-inputfile=tmp_input.txt -outputfile=output.txt -subtree=yes"; # > $log_file";
  
  my $nbTrees = 0 ; # scalar @trees_tab - 1;
  if($bootstrap eq "yes" || $consensus eq "yes" ){
  #if($bootstrap eq "yes"){
	#$nbTrees -= 1; 
  }
  
  #== The program need at least 2 trees
  if((scalar @trees_tab < 2)){
    exit_program(-1,$return_file,"PERL : nombre d'arbres invalide");
  }
  
  print_minidoc(); 
  
   # my $nbHGT=0;
   # my $nbHGT2=0;
	
  for (my $i=0;(($i< scalar @trees_tab) && $trees_tab[$i] =~ ";");)
  {
    #print STDOUT "\n==== $i\n";
    open(IN,">$tmp_input") ||  die "Cannot open $tmp_input";
	
    #================================================================================
	#= In the bootstrap case, we need to change the speciesroot and generoot option
	#= for "file" from the first replicate.
	#================================================================================
	#if($bootstrap eq "yes"){
	if($bootstrap eq "yes" || $consensus eq "yes" ){
		if($i == 0){ 
			print IN $trees_tab[0] . $trees_tab[1];
			$i=2;
			$cmd .= " -randbootstrap=$rand_bootstrap";
        }
		else{
			print IN $trees_tab[0] . $trees_tab[$i++];
		}
		#print STDOUT $trees_tab[0] . $trees_tab[$i];
  		if($i > 2){
			if($cmd !~ /printWeb=no/){
				$cmd .= " -printWeb=no ";
			}
			
		}
		if($nbTrees == 1){
			$cmd =~ s/-generoot=[a-z][a-z]* //;
			$cmd =~ s/-speciesroot=[a-z][a-z]* //;
			if($consensus eq "yes"){
				$cmd .= " -generoot=midpoint -speciesroot=midpoint";
			}
			else{
				$cmd .= " -generoot=file -speciesroot=file";
			}
			#print STDOUT "PERL : 1 : nbTress=$nbTrees";
			
			#= linux
			`cp $generootfile $generootfiletmp`;
			`cp $speciesrootfile $speciesrootfiletmp`;			
			
			#= windows
			# `copy $generootfile $generootfiletmp`;
			# `copy $speciesrootfile $speciesrootfiletmp`;
		}
		if($nbTrees > 1){
			#== linux like
			`cp $generootfiletmp $generootfile`;
			`cp $speciesrootfiletmp $speciesrootfile`;
			
			#== windows
			# `copy $generootfiletmp $generootfile`;
			# `copy $speciesrootfiletmp $speciesrootfile`;
		}
		
		
		
		if($rand_bootstrap == 1){
			$rand_bootstrap = 0;
			$cmd =~ s/randbootstrap=1/randbootstrap=0/;
		}
		else{
			$rand_bootstrap = 1;
			$cmd =~ s/randbootstrap=0/randbootstrap=1/;
		}
    }
    else{
        print IN $trees_tab[$i++] . $trees_tab[$i++];
    }
    close(IN);
    
	#== traitement des cas de transferts avant detection
	if($stepbystep eq "yes"){
		if(-e $outputWeb){
			my @prehgt;
			my @prehgt_valide;
			my @prehgt_sequence;
			open(PREHGT,"$outputWeb") || die "Probleme avec l'ouverture du fichier $outputWeb";
			my @tab_prehgt = <PREHGT>;
			foreach my $ligne( @tab_prehgt){
				chomp($ligne);
				if($ligne =~ /^hgt_reels =/){
					$ligne =~ s/ = /,/;
					@prehgt = split(",",$ligne);
				}
				if($ligne =~ /^hgt_valide =/){
					$ligne =~ s/ = /,/;
					@prehgt_valide = split(",",$ligne);
				}
				if($ligne =~ /^hgt_sequence/){
					$ligne =~ s/ = /,/;
					@prehgt_sequence = split(",",$ligne);
				}
			}
			close(PREHGT);
			print STDOUT join("-",@prehgt) . "\n";
			print STDOUT join("-",@prehgt_valide) . "\n";
			print STDOUT join("-",@prehgt_sequence) . "\n";
			
			open(PREHGT,">$prehgtfile") || die "Impossible d'ouvrir $prehgtfile";
			print STDOUT "\nnbHGT  = " . $#prehgt_valide . "\n\n";
			
			#exit;
			
			for(my $i=1;$i< scalar(@prehgt_valide);$i++){
				print PREHGT $prehgt_sequence[$i] . " " . $prehgt[4*$i-3] . " " . $prehgt[4*$i-2] . " " . $prehgt[4*$i-1] . " " . $prehgt[4*$i-0] . " " . $prehgt_valide[$i] . "\n";
			}
			close(PREHGT);
		}
	}
	
#exit;
	
	print STDOUT "\nComputation " . ++$nbTrees . " in progress...";
	
	if($consensus eq "yes"){
		#$cmd .= " -bootstrap=yes";  
		$cmd =~ s/bootstrap=no/bootstrap=yes/;
	}

	
	print STDOUT "\nCommande=$cmd";
    execute_hgt("$cmd >> $log_file");
	if ($consensus eq "yes"){
		open(INFO, $outputWeb);	# Open the file
		my @lines = <INFO>;		# Read it into an array	
		close(INFO);			# Close the file
		foreach my $line( @lines){
			
			if ($line =~ /hgt = /){
				if ($premier ==1){
					$premier = 0;				
				}
				else{
					$line=~ s/hgt = /,/;
				}
				$line =~ s/\s+$//;   
				$hgt = $hgt . $line;				
			}
		}
	}
	
	if( ! -e $results){
	  print STDOUT "\n\nRUN_HGT : An error has occured during computation. Check the log file ($log_file) for more details ! ";
	  exit -1;
    }
	print STDOUT "\nformatting results...";
	
    if((($i == 2) && ($viewtree eq "yes")) || (($i == 1) && ($viewtree eq "yes"))){
       exit_program(0,$return_file,"RUN_HGT : We just want to see the input tress");
    }
    
  #===========================================================================
  #========================= LECTURE DES RESULTATS ===========================
  #===========================================================================
    open(IN,"$results") ||  die "\nCannot open $results ! !!";
    @hgt_tab = <IN>;
    close(IN); 
 	
    exit_program(-1,$return_file,"RUN_HGT : result file empty") if(scalar @hgt_tab == 0);  
     
    $mode = $hgt_tab[1];
    chomp($mode);
   
    if(  ($i == 2) || (($bootstrap eq "no") && ($consensus eq "no"))  ){
		
		print STDOUT "done";
        $nbHGT=0;
        $nbHGT2=0;
		
        for(my $j=2;$j<= scalar @hgt_tab;){
			if($hgt_tab[$j] =~ /^[0-9]/){
                my $cpt = read_line($hgt_tab[$j++]); #2
				
				for(my $k=0;$k<$cpt;$k++){ 
                    my $hgt_number = read_line($hgt_tab[$j++]); # HGT 1 / 2 
					
					# my $hgt_number = read_line($hgt_tab[1]);
					# my $hgt_number = read_line($hgt_tab[2]);
					
                    my $source_list = read_line($hgt_tab[$j++]); # D
                    my $dest_list = read_line($hgt_tab[$j++]); # B
                    
                    $source_list = join(", ",sort split(", ",$source_list));
                    $dest_list = join(", ",sort split(", ",$dest_list));
                    
					my $transfer_description2 = read_line($hgt_tab[$j++]); # From branch D--6 to branch B--7
					my $transfer_description = "From subtree ($source_list) to subtree ($dest_list)"; 
					
                    my $criterion_list = read_line($hgt_tab[$j++]); # RF = 2 , LS = 0.444 , BD = 1.0
					# if($mode eq 'mode=multicheck'){
						my $criterion_list2 = read_line($hgt_tab[$j++]); # rRF = 2 , rLS = 0.444 , rBD = 1.0
					# }
				 
                    $hgt_number_tab{"$source_list -> $dest_list"}      = $hgt_number;
                    $hgt_description_tab{"$source_list -> $dest_list"} = $transfer_description;
                    $hgt_compteur_tab{"$source_list -> $dest_list"}    = 1;
                    $hgt_compteur_tab_reverse{"$dest_list -> $source_list"}    = 1;
                    $hgt_criterion_tab{"$source_list -> $dest_list"}   = $criterion_list;
					$hgt_criterion_tab2{"$source_list -> $dest_list"}   = $criterion_list2;
                    $hgt_nbHGT_tab{"$source_list -> $dest_list"}       = $cpt;
                    $hgt_pos[$nbHGT++] = "$source_list -> $dest_list";
                }
                # if($mode eq 'mode=multicheck'){
                    $hgt_pos2[$nbHGT2++] = read_line($hgt_tab[$j++]); #read_line($hgt_tab[$j++]);
                    
                # }
				
				
            }
            else{
				
				@tmp_tab_init = split(",",$hgt_tab[0]);
                @tmp_tab = split(" ",$hgt_tab[$j]);
                $total_hgt = $tmp_tab[1];
                $total_trivial = $tmp_tab[2];
				
				
                if( ($bootstrap eq "no") && ($consensus eq "no") ){
					
    				open(OUT,">>$outputfile") || die "Cannot open $outputfile";
                    if($total_hgt > 0){
                        print_result();
                    }
                    else{
                        print OUT " : no HGTs have been found !";
                    }
                   
                    close(OUT);
                    
                    exit_program($val_retour,$return_file,"PERL : pas de bootstrap, on traite un seul input");
                    %hgt_number_tab=();
                    %hgt_description_tab=();
                    %hgt_compteur_tab=();
                    %hgt_criterion_tab=();
					%hgt_criterion_tab2=();
                    %hgt_nbHGT_tab=();
                    @hgt_pos=();
                    @hgt_pos2=();
                } 
                $j = (scalar @hgt_tab) + 1;
            }
        } 
		if ($bootstrap eq "yes")
		{
			copy($results, $firstResults) or die "File cannot be copied.";
			copy($outputWeb, $firstOutputWeb) or die "File cannot be copied.";
			copy($results2, $firstResults2) or die "File cannot be copied.";
		}
    }


    if(  (($bootstrap eq "yes") || ($consensus eq "yes"))   && ($i > 2))
	{
		print STDOUT "done";
		if($bootstrap eq "yes") {
			$nbHGT=0;
			$nbHGT2=0;
		}
        for(my $j=2;$j<= scalar @hgt_tab;){
            if($hgt_tab[$j] =~ /^[0-9]/){
                my $cpt = read_line($hgt_tab[$j++]);
                for(my $k=0;$k<$cpt;$k++){ 
					my $hgt_number = read_line($hgt_tab[$j++]);
                    my $source_list = read_line($hgt_tab[$j++]);
                    my $dest_list = read_line($hgt_tab[$j++]);
                    $source_list = join(", ",sort split(", ",$source_list));
                    $dest_list = join(", ",sort split(", ",$dest_list));
					my $transfer_description2 = read_line($hgt_tab[$j++]);
					my $transfer_description = "From subtree ($source_list) to subtree ($dest_list)"; 
					
                    my $criterion_list = read_line($hgt_tab[$j++]);
                    my $criterion_list2 = read_line($hgt_tab[$j++]);
					
                    if(exists $hgt_compteur_tab{"$source_list -> $dest_list"}){
                        $hgt_compteur_tab{"$source_list -> $dest_list"}   += 1;
						# print STDOUT "\nJe rentre dans le if exist " . $transfer_description . $hgt_compteur_tab{"$source_list -> $dest_list"} ;						
                    }
					else{
						if($consensus eq "yes"){
								$hgt_number_tab{"$source_list -> $dest_list"}      = $hgt_number;
								$hgt_description_tab{"$source_list -> $dest_list"} = $transfer_description;
								$hgt_compteur_tab{"$source_list -> $dest_list"}    = 1;
								$hgt_criterion_tab{"$source_list -> $dest_list"}   = $criterion_list;
								$hgt_criterion_tab2{"$source_list -> $dest_list"}   = $criterion_list2;
								$hgt_nbHGT_tab{"$source_list -> $dest_list"}       = $cpt;
								$hgt_pos[$nbHGT++] = "$source_list -> $dest_list";
								# print STDOUT "\nJe rentre dans le if consensus " . $transfer_description . $hgt_compteur_tab{$hgt_pos[$nbHGT-1]} ;
								# print STDOUT "\nJe rentre dans le if consensus " . $transfer_description . $hgt_compteur_tab{"$source_list -> $dest_list"} ;
						}	
						if(exists $hgt_compteur_tab_reverse{"$source_list -> $dest_list"}){
                       	 	$hgt_compteur_tab_reverse{"$source_list -> $dest_list"}   += 1;
							#print STDOUT "\nJe rentre dans le if exist " . $transfer_description . $hgt_compteur_tab_reverse{"$source_list -> $dest_list"} ;						
                    	}
					}					
                }
				# if(($mode eq 'mode=multicheck') && ($consensus eq "yes")){
				if($consensus eq "yes"){
					  $hgt_pos2[$nbHGT2++] = read_line($hgt_tab[$j++]);
				}	
                # elsif($mode eq 'mode=multicheck'){
				else{
                    $j++;
                }
            }
            else{		
                $j = (scalar @hgt_tab) + 1;			
            }
        } 
		if($consensus eq "yes")
		{
			$total_hgt =  $nbHGT;
		}
		
		print STDOUT "\nLe nombre de transferts" . $nbHGT ;
    }
} 

if ($bootstrap eq "yes")
{
	# On recopie les premiers résultats sauvegardés.
	copy($firstResults, $results) or die "File cannot be copied.";
	copy($firstOutputWeb, $outputWeb) or die "File cannot be copied.";
	copy($firstResults2, $results2) or die "File cannot be copied.";
	unlink($firstResults);
	unlink($firstOutputWeb);
	unlink($firstResults2);
}

if($consensus eq "yes"){
	
	open(OUT,">>$outputfile" ) || die "Cannot open $outputfile";
	print_result_consensus();
	close(OUT);
}

if($bootstrap eq "yes"){
	open(OUT,">>$outputfile") || die "Cannot open $outputfile";
	print_result();
	close(OUT);
}

exit_program($val_retour,$return_file,"PERL : fin normale du programme");
              
#===============================================================================
#=============================== FUNCTIONS =====================================
#===============================================================================

sub read_line{
  my ($line) = @_;
  chomp($line);
  return $line;
}

sub exit_program{
  my($val,$file,$message) = @_;
  open(RET,">$file") || die "Cannot open $file";
  print RET $val;
  close(RET);
  print STDOUT "\n";
  #print STDOUT "\nexit=>$message";
  &clean();
  exit;
}

sub execute_hgt{
    my ($cmd) = @_;
    my $retour = 0;
    $retour= system($cmd); #print STDOUT system($cmd);	
}


sub print_result{
    my $sim_hgt = "yes";
    open (SIM , ">res_sim.txt") if ($sim_hgt eq "yes");
    my $nbHGT2=0;
    my $newGroup=0;
    my $cpt=1;
    my @tmp_tab = split(",",$hgt_tab[0]);
	
	print OUT	"==================================================================================\n";
	print OUT	"| Program : HGT Detection 3.4 - March, 2012                                      |\n";
    print OUT   "| Authors   : Alix Boc, Alpha Boubabcar Diallo and Vladimir Makarenkov           |\n";
    print OUT   "|             (Universite du Quebec a Montreal)                                  |\n";
	print OUT	"| This program computes a unique scenario of horizontal gene transfers (HGT) for |\n"; 
    print OUT   "| the given pair of species and gene phylogenetic trees.                         |\n";
	print OUT	"==================================================================================\n";
	
	
    print OUT "\nSpecies tree :\n". $trees_tab[0] . "\nGene Tree :\n" . $trees_tab[1];
                   
    print OUT "\n\n=============================================";
	  print OUT "\n= Criteria values before the computation ";
	  print OUT "\n=============================================";	
	if($bootstrap eq "yes"){
	  printf (OUT "\nRobinson and Foulds distance (RF) = %d",$tmp_tab_init[0]);
	  printf (OUT "\nLeast-squares coefficient(LS)     = %1.3lf",$tmp_tab_init[1]);
	  printf (OUT "\nBipartition dissimilarity         = %1.1lf\n",$tmp_tab_init[2]);
	}
	else{
	  printf (OUT "\nRobinson and Foulds distance (RF) = %d",$tmp_tab[0]);
	  printf (OUT "\nLeast-squares coefficient(LS)     = %1.3lf",$tmp_tab[1]);
	  printf (OUT "\nBipartition dissimilarity         = %1.1lf\n",$tmp_tab[2]);
	}
	
	printf(OUT "\n\nBootstrap values were computed with %d gene trees",$nbTrees) if($bootstrap eq "yes");

 #   print STDOUT "mode=$mode";
 
    print OUT "\n\n";
    foreach my $elt( @hgt_pos){
		 
        # if(($newGroup == 0 ) && ($mode eq 'mode=multicheck')){
		if($newGroup == 0 ) {
            print OUT "\n================================================================"; 
            if($hgt_nbHGT_tab{"$elt"} == 1){
                print OUT "\n| Iteration #$cpt : ". $hgt_nbHGT_tab{"$elt"} ." HGT was found";
                     
                #print OUT "\n| [" . $hgt_nbHGT_tab{"$elt"} . " HGT was found in this iteration]";
            }
            else{
              print OUT "\n| Iteration #$cpt : ". $hgt_nbHGT_tab{"$elt"} ." HGTs were found";
               # print OUT "\n| [" . $hgt_nbHGT_tab{"$elt"} . " HGTs were found in this iteration]";
            }
            
            print OUT "\n================================================================";
            print OUT "\n|";
            $newGroup = 1;
            $cpt++; 
        }
        
        print OUT "\n| "  . $hgt_number_tab{"$elt"};
		
		my @tmp = split(" -> ",$elt);
		my $elt_reverse = $tmp[1] . " -> " . $tmp[0];
		
		#print OUT "\n$elt\n$elt_reverse";
		my $hgt_type = "Trivial";
		if(($bootstrap eq "yes")) { #&&($hgt_number_tab{"$elt"} !~ "Trivial")){
			if($hgt_number_tab{"$elt"} !~ "Trivial"){
				print OUT " Regular";
				$hgt_type = "Regular";
			}
			printf(OUT " (bootstrap value = %3.1lf%% inverse = %3.1lf%%) ",$hgt_compteur_tab{$elt}*100/$nbTrees,$hgt_compteur_tab_reverse{$elt_reverse}*100/$nbTrees); 
			my @tmp_sim_tab = split("->",$elt);
			printf(SIM "%s<>%s<>%3.1lf<>%3.1lf<>%s\n",$tmp_sim_tab[0],$tmp_sim_tab[1],$hgt_compteur_tab{$elt}*100/$nbTrees,$hgt_compteur_tab_reverse{$elt_reverse}*100/$nbTrees,$hgt_type) if ($sim_hgt eq "yes"); 
		}
		print OUT "\n| "  . $hgt_description_tab{"$elt"};
        print OUT "\n| "  . $hgt_criterion_tab{"$elt"};
		print OUT "\n| "  . $hgt_criterion_tab2{"$elt"};
        if($mode eq 'mode=monocheck'){
            print OUT "\n================================================================\n"; 
        }
        else{
            print OUT "\n| "; 
        }
        my $tmp = "HGT " . $hgt_nbHGT_tab{"$elt"} . " / " . $hgt_nbHGT_tab{"$elt"} . " ";
        my $tmp2 = $tmp . " Trivial ";
        
       if(( $tmp =~ $hgt_number_tab{"$elt"}) ||($tmp2 =~ $hgt_number_tab{"$elt"})){
            if ($mode eq 'mode=multicheck') {
				print OUT "\n================================================================";
				print OUT "\n| After this iteration the criteria values are as follows :";
				print OUT "\n| " . $hgt_pos2[$nbHGT2++] ;
				print OUT "\n================================================================\n";
			}
				$newGroup=0; 
        }
    } 
    
    print OUT "\nTotal number of HGTs : $total_hgt ";
    print OUT "(". ($total_hgt-$total_trivial) ." regular + " . $total_trivial . " trivial HGTs)" if( $total_trivial > 0);
    
    $val_retour = $total_hgt;
	
	open(OUTWEB,">>$outputWeb");
	if($bootstrap eq "yes"){
		print OUTWEB "\nbootHGT = ";
		my $first=0;
		foreach my $elt( @hgt_pos){
			
			if(($bootstrap eq "yes")&&($hgt_number_tab{"$elt"} !~ "Trivial")){
				if($first==0){
					printf(OUTWEB "%3.0lf",$hgt_compteur_tab{$elt}*100/$nbTrees);
				}
				else{
					printf(OUTWEB ",%3.0lf",$hgt_compteur_tab{$elt}*100/$nbTrees);
				}
				$first=1;
			}
			elsif($bootstrap eq "yes"){
			    if($first==0){
					printf OUTWEB "0";
				}
				else{
					printf OUTWEB ",0";
				}
				$first=1;		
			}
		}
	}
	close OUTWEB;
	
    close (SIM) if ($sim_hgt eq "yes");
}



sub par_elt { return $hgt_compteur_tab{"$b"}  <=> $hgt_compteur_tab{"$a"} }

sub print_result_consensus{
    my $nbHGT2=0;
    my $newGroup=0;
    my $cpt=0;
    my @tmp_tab = split(",",$hgt_tab[0]);
    
	my @out =  sort par_elt @hgt_pos;
	
	print OUT	"==================================================================================\n";
	print OUT	"| Program : HGT Detection 3.4 - March, 2012                                      |\n";
    print OUT   "| Authors   : Alix Boc, Alpha Boubabcar Diallo and Vladimir Makarenkov           |\n";
    print OUT   "|             (Universite du Quebec a Montreal)                                  |\n";
	print OUT	"| This program computes a unique scenario of horizontal gene transfers (HGT) for |\n"; 
    print OUT   "| the given pair of species and gene phylogenetic trees.                         |\n";
	print OUT	"==================================================================================\n";
	

	open(OUTCONS,">$outputConsensus") || die ("Erreur d'ouverture de outputConsensus") ;
	
	foreach my $elt( @out){
		printf(OUT "\n(consensus value = %d  of %d = %3.0lf%) ",$hgt_compteur_tab{"$elt"}, $nbTrees,  (($hgt_compteur_tab{"$elt"} * 100 ) / $nbTrees )); 
		print OUT "\n "  . $hgt_description_tab{"$elt"};
		if ( (($hgt_compteur_tab{"$elt"} * 100 ) / $nbTrees )>= $cons){		
			# (my @tab_boot) = ($chaine =~ /\)([0-9][0-9]*):/g);
			(my @tab_cons) = ($hgt_description_tab{"$elt"} =~ /\([^)]+\)/g);
			foreach my $espece( @tab_cons){	
					$nb = ($espece =~ tr/,//) + 1 ;
					$espece =~ s/\(/ /;
					$espece =~ s/\)/ /;
					$espece =~ s/,/ /g;				
					printf OUTCONS "\n". $nb . $espece;
			}		
		}
    } 
	

	
	my $cmd_cons = "hgt_find_branches -inputfile=$inputfile_bkp -path=$path";  
	$cmd_cons          = "usagers/" . $cmd_cons;
	execute_hgt("$cmd_cons >> $log_file1");
    close OUTCONS;
	$val_retour = scalar(@out);
    print OUT "\n\n\nTotal number of HGTs : $val_retour ";
    # print OUT "(". ($total_hgt-$total_trivial) ." regular + " . $total_trivial . " trivial HGTs)" if( $total_trivial > 0);
    
    
	# $val_retour = 20;
	open(OUTWEB,">>$outputWeb");
	printf(OUTWEB "\nbootHGT=");
	my $first=0;
	foreach my $elt( @out){
		if($first==0){
			printf(OUTWEB "%3.0lf",$hgt_compteur_tab{$elt}*100/$nbTrees);
		}
		else{
			printf(OUTWEB ",%3.0lf",$hgt_compteur_tab{$elt}*100/$nbTrees);
		}
		$first=1;
    } 
	close OUTWEB;
}

sub print_title{
	print STDOUT "===================================================================================================\n";
	print STDOUT "| HGT-DETECTION V.3.4 (March, 2012) by Alix Boc, Alpha Boubacar Diallo and Vladimir Makarenkov |\n"; 
	print STDOUT "===================================================================================================\n";
}

sub print_minidoc{
	print STDOUT "\nCheck the file $log_file for the computation details";
	print STDOUT "\nCheck the file $outputfile for the program output\n";
}

sub print_description{
	print STDOUT	"==========================================================================================================\n";
	print STDOUT	"| Program : HGT Detection 3.4 - March, 2012                                                           |\n";
    print STDOUT    "| Authors   : Alix Boc, Alpha Boubacar Diallo and Vladimir Makarenkov (Universite du Quebec a Montreal)  |\n";
	print STDOUT	"| This program computes a unique scenario of horizontal gene transfers (HGT) for                         |\n"; 
    print STDOUT    "| the given pair of species and gene phylogenetic trees.                                                 |\n";
	print STDOUT	"==========================================================================================================\n";
}

sub print_help{
	print STDOUT "\nUsage :\nperl run_hgt.pl -inputfile=[inputfilename] -outputfile=[outputfilename] -criterion=[rf|ls|bd]";
	print STDOUT "-speciesroot=[midpoint|file] -generoot=[midpoint|file|bestbipartition]";
	print STDOUT "-scenario=[unique|multiple] -nbhgt=[maxhgt] -path=[path] -bootstrap=[no|yes]";
	print STDOUT "\n\nsee README.txt file for more detail.";

}

sub clean(){
	`rm -rf errorFile.txt *~ inputfileFormated.txt input_.txt tmp_input.txt geneRoot* speciesRoot* nomorehgt.txt res_sim.txt`;
}
