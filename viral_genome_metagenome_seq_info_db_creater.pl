#!/usr/local/env perl
#####################################################################################################################
#This program populates a database that has already been desigened and created by another script.  The program      #
#knows about the database, it accessses the database by asking permitted user's login and password.  The script     #
#knows the sorce of data, does data abstraction, mining, processing, modification, transporation and laods the      #
#crossed checked data to exact tables and location in the database.  The program constantly verifies that the       #
#correct and right data is entering the database at in the right table and in right column.  The program extracts   #
# gene and contif sequences from different public data base, cross checks the validity of protein sequence in       #
#6 frames, and the translated sequences from ORFs are validated with the NCBI's protein database before transporting#
#data into the database.                                                                                            #
#####################################################################################################################

use lib '/usr/local/projects/GOSII/skulshre/bioperl-1.4';
use Bio::Perl;
use Bio::Tools::CodonTable;
use lib '/usr/local/projects/GOSII/skulshre/';
use BeginPerlBioinfo;

##############################################################################################################


# for DB
use lib '/usr/local/projects/GOSII/skulshre';
use DBI;
use DBD::mysql;

my $sqlhost = "mysql.tigr.org";
my $database = "viraldb";
my $user = "skulshre";
my $password = "10355Science!";


print "connecting to $sqlhost...";	
my $connect_to_host = DBI->connect("DBI:mysql:$database:$sqlhost", $user, $password) or die "Can't connect to MySql server at $sqlhost\n";
print "succeeded.\n";

################################################################################################################


#THE NAME OF THIS FILE IS:/home/skulshre/sk_mylocal_programs/sk_files/genbank_all_virus_genome_seq_parsed.pl


# READING DIRECTORY FOR ALL THE FILES WITH FULL PATH FOR EACH FILE IN ONE ARRAY
#FILE IN THIS DIRECTORY testfile1.fa HAS BEEN TESTED ON THIS PARSER FILE: sdsu_phage_genome_seq_parsed.pl
#THIS FILE LOCATION HAS BEEN CHANGED TO FOLLOWING LOCATION ON 06/15/08
#@phagefiles = </home/skulshre/sk_mylocal_programs/sk_files/test_dir/*>;

#A NEW DIRECTORY LOCATION TO RUN PROGRAMS AND FILES HAS BEEN CHANGED ON 06/16/2008 TO:
@phagefiles = </usr/local/projects/GOSII/skulshre/sk_mylocal_programs/sk_files/test_dir/*>;


$genome_count = 0;
$genomecount2 = 0;
$contig_seq;
$percent_gc;
%contig_all_features_hash=();
$total_db_cds_count=0;
$total_db_cdsseqs_ct2 = 0;
$db_protein_count  = 0;
$total_db_protseqs_ct2 = 0;


foreach my $name (@phagefiles ){
    
    next if ($name =~ /\~$|\#$|^\#/);
    
    if ($name =~ /\S+.*(\.gbff$)/ ){ 

# FOR VIRAL DB ........................................................................
	
	($genomecount2, $total_db_cdsseqs_ct2,$total_db_protseqs_ct2) = contigSeq($name);
	
# }
	
	print "global genome count value = ", $genome_count, "\n", "genome count value from sub = ",$genomecount2, "\n" ;
	print "global total db CDS seqs count = ", $total_db_cds_count, "\n", "total cds count value returned from sub = ",$total_db_cdsseqs_ct2, "\n" ;
	
	print "global Total db proteins count = ", $db_protein_count, "\n", "total protein count returned from sub = ", $total_db_protseqs_ct2, "\n";
	
    }elsif ($name =~ /\S+.*(metagenomicsdb$)/){

# FOR METAGENOMICS DB ....................................................................
	
	
#	system("perl /usr/local/projects/GOSII/skulshre/sk_mylocal_programs/sk_files/viral_metagene3.pl > /usr/local/projects/GOSII/skulshre/out1 &");
	
%hash_taxonomy = (BEARPAW_SMPL_20031007 => 'Viruses; Hotspring; Yellowstone; Bearpaw',
	       GIS_HUMAN_FECAL_VIRUS_RNA_SMPL => 'Viruses; Fecal; Human; RNA',
	       OCTOPUS_SMPL_20031004 => 'Viruses; Hotspring; Yellowstone; Octopus',
	       SDSU_HUMAN_BLOOD_VIRUS_SMPL => 'Viruses; Blood; Human; SDSU',
	       SDSU_HUMAN_FECAL_VIRUS_SMPL => 'Viruses; Fecal; Human; DNA',
	       SDSU_MARINE_SEDIMENT_VIRUS_SMPL => 'Viruses; Marine; Sediment; SDSU',
	       SDSU_MISSION_BAY_VIRUS_SMPL => 'Viruses; Marine; Planktonic; Mission bay',
	       SDSU_SCRIPPS_PIER_MARINE_VIRUS_SMPL => 'Viruses; Marine; Planktonic; Scripps pier',
	       UBC_MARINE_VIRAL_RNA_SMPL => 'Viruses; Marine; Planktonic; Marine viral RNA',
	       UDEL_CHESAPEAKE_VIROPLANKTON_SMPL => 'Viruses; Marine; Planktonic; Chesapeak',
	       ULE_EQUINE_FAECAL_PHAGE_SMPL => 'Viruses; Fecal; Equine; Equine fecal phage'
    );

%hash_taxon_id = (Virus => '5',
		  Hotspring => '1',
		  Fecal => '2',
		  Blood => '3',
		  Marine => '4',
		  Yellowstone => '1',
		  Human => '1',
		  Equine => '2',
		  Planktonic => '1',
		  Sediment => '2',
    );

my $taxonomy;
my $form = 'metagenome';

$main_dir =  '/usr/local/projects/GOSII/external_metagenomic_projects';

@arrray_metagene_dir = </usr/local/projects/GOSII/external_metagenomic_projects/*>;
print "array length =", $length = @arrray_metagene_dir,"\n";

print "directories are :\n@arrray_metagene_dir\n";


foreach $my_metagenedir_name (@arrray_metagene_dir ){

    $my_metagenedir_name =~ s/external_metagenomic_projects/skulshre\/external_metagenomic_projects/;
    
    if (-d $my_metagenedir_name  ){
	print "I found the directory : $my_metagenedir_name\n";
	($make_dir = $my_metagenedir_name) =~ s/^\S.*projects(\/\S+.*)$/$1/;
	$main_dir .= $make_dir;

#	system ("cp -rf $main_dir/clr_range_filter_orf $my_metagenedir_name/");
#	system ("cp -rf $main_dir/clr_range_filter_pep $my_metagenedir_name/");
	
	$main_dir = '/usr/local/projects/GOSII/external_metagenomic_projects';

	$make_dir =~ s/\///g;
	my $contig_name = $make_dir;

	next;
	
    }else{
	print "I am making directory :", $my_metagenedir_name, "\n";
	system("mkdir $my_metagenedir_name");
	system ("cp -rf $main_dir/clr_range_filter_orf $my_metagenedir_name/");
	system ("cp -rf $main_dir/clr_range_filter_pep $my_metagenedir_name/");

	print " I made the directory: $my_metagenedir_name ", "\n"; 

	$main_dir = '/usr/local/projects/GOSII/external_metagenomic_projects';
	$make_dir =~ s/\///g;
	my $contig_name = $make_dir;
	
    }
}


# /usr/local/projects/GOSII/skulshre/external_metagenomic_projects/BEARPAW_SMPL_20031007/clr_range_filter_orf/1113077373523_default 

my $orflength = 0;
my $metagene_contig = 0;
my $avg_size = 0;
my $orf='';
my $taxon_id='';
my $lnth='';
my $orfcount=0;
my @contigcdsseqs=();
my $pepseq='';
my $protein_count=0;
my $protein_name='';
my $pep_signal = 0;
my $pep_begin_signal=0;
my $contig_total_proteins=0;

my @metagene_contigs = </usr/local/projects/GOSII/skulshre/external_metagenomic_projects/*>;  

foreach my $metagene_contig_name (@metagene_contigs){ # 1
    
    print "the metagene contig : $metagene_contig_name\n";
    ($contig_name = $metagene_contig_name) =~ s/^\/usr\S.*\/(\S.*$)/$1/;
    
    $species = $contig_name;
    $taxonomy = $hash_taxonomy{$contig_name};
    $taxonomy_copy = $taxonomy;
    $taxonomy_copy =~ s/\s//g;

    ($a, $b, $c, $d ) = split(/\;/, $taxonomy_copy);
    $taxon_id = $hash_taxon_id{$a}; 
    $taxon_id .= $hash_taxon_id{$b}; 
    $taxon_id .= $hash_taxon_id{$c}; 

 #   $taxon_id = $taxon_id_a.$taxon_id_b.$taxon_id_c;

    ($read_contig_file )= $metagene_contig_name =~ /metagenomic_projects\/(\S+.*$)/;
    
    print "path is : $read_contig_file\n";
    
    @subdir = </usr/local/projects/GOSII/skulshre/external_metagenomic_projects/$read_contig_file/*>;
    
# GOING INTO EACH METAGENE DIRECTORY
    foreach my $orfdir (@subdir){  # 2
#	next if ($orfdir =~ /(pep)$/);
	


# GOING INTO ORF DIRECTORY
	if ($orfdir =~ /(orf)$/){ # 3
	    print "1= dir is: ", $orfdir, "\n";
	    @orf_dir = <$orfdir/*>;
	    
# GOINT INTO SUB DIR
	    foreach my $orf_subdir (@orf_dir ) { # 4
		print "2= the next level dir is :", $orf_subdir, "\n";
		
		@orf_files = <$orf_subdir/*>;
		
		print "3= @orf_files\n";
		
# ACCESSING .fna FILE
		foreach my $orf_file (@orf_files ){# 5
		    next if ($orf_file !~ /\S+.*(\.fna)$/);
		    print "4= orf file is :", $orf_file, "\n";
		    
# READING .fna FILE		   
		    if ($orf_file =~ /\S+.*(\.fna)$/){ # 6
		    open (READFILE, "<$orf_file") or die "could not open orf file : $orf_file\n" ;
		    $orfcount = 0;
		    while ($line = <READFILE>){
			
			if(($line =~ /^([A-Z]+)$/) && ($line !~ /^\>/ )){
			    $orf .= $line;
			    $orfsignal = 1;
			   

			}elsif(($line =~ /^\>/) && ($orfsignal == 1 )){
			    $orfcount += 1;
			    $orf =~ s/''|\W|\n|\s//g;
			    $lnth = length $orf;
			    push (@contigcdsseqs, $lnth.":::".$orf."***".$orfcount);
			    print "the orf after push is:\n", $lnth.":::".$orf."***".$orfcount, "\n";
			    $orfsignal = 0;
			    $metagene_contig .= $orf;
			    $orf='';

			}
			
		    } # while loop
		    
# PROCESSING LAST ORF OUTSIDE THE LOOP
		    $orf =~ s/''|\W|\n|\s//g;
		    push (@contigcdsseqs, $lnth.":::".$orf."***".$orfcount);
		    $orfsignal = 0;
		    $metagene_contig .= $orf;
		    $orf='';
		    $lth = @contigcdsseqs;
		    $orfcount += 1;

		    print "the size of the cds array is =", $lth, "\n";
		    print "the cds count after loop = ", $orfcount, "\n";

		    $metagene_contig =~ s/''|\W|\n|\s//g;
		    $contig_gc= ($metagene_contig  =~ tr/GCgc//);
		    $contig_gc = ($contig_gc/length $metagene_contig)*100;
		    
		    foreach my $cdsinfo (@contigcdsseqs ){
			
			($size, $cdsseq_cdscount) = split(/:::/, $cdsinfo);
			
			next if ($size =~ /''|NULL/);
			$avg_size += $size; 	
		
		    }
		    
		    $contig_size = $avg_size/($length =@contigcdsseqs );
		    $avg_size=0;
		    creatMetcontigTable($contig_name, $species, $taxonomy, $form, $contig_gc, $contig_size, $taxon_id);   

		    }# 6 if .fna file found
		} # 5
	    } # 4
	    print "\n";
	    
	}# 3
	print "protein count after contigs file is read = ", $protein_count, "\n";
	$protein_count = 0;

# GOING INTO PEP SUBDIR
	if ($orfdir =~ /(pep)$/){# 7

	    $pepdir = $orfdir;

	    print "5= dir is: ", $pepdir, "\n";
	    @pep_dir = <$pepdir/*>;

	    foreach my $pep_subdir (@pep_dir ) { # 8
		print "6= the next level dir is :", $pep_subdir, "\n";		
		@pep_files = <$pep_subdir/*>;
		print "7= @pep_files\n";

# ACCESSING .faa FILE
		foreach my $pep_file (@pep_files ){# 9
		    next if ($pep_file !~ /\S+.*(\.faa)$/);
		    print "8= pep file is :", $pep_file, "\n";
		    
		    if ($pep_file =~ /\S+.*(\.faa)$/){ # 10
			open (READFILE, "<$pep_file") or die "could not open pep file : $pep_file\n" ;
			
			$pep_begin_signal = 0;
			while ($line = <READFILE>){ # 11
			    
			    if(($line =~ /^\>/)&& ($pep_begin_signal == 0)){
			       
				print "\npep signal *= ", $pep_signal, "\n", $line, "\n";
				@prot_names = split (/\s/,$line);
				$protein_name = @prot_names[0];
				$protein_name =~ s/\>//g;
				print "protein name for first protein = ", $protein_name, "\n",
				#$pep_signal =  1;
				$pep_begin_signal = 1;
				$protein_count+= 1;
				next;
			    }
			    elsif(($line =~ /^\>/)&& ($pep_signal == 1)){
			
				print "values for protein table are:\n",
				"protein name = ", $protein_name, "\n",
				"contig name = ", $contig_name, "\n",
				"protein sequence:\n",
				$pepseq."***".$protein_count, "\n";
				print "protein count before protein table = ", $protein_count, "\n";

				my $protein_id = creatMetProteinTable($pepseq, $contig_name, $protein_name);

				$cds = @contigcdsseqs[$protein_count-1];
				
				($cdssize, $cdsseq_cdscount)  = split(/:::/, $cds);
				($cdssequence, $count) = split(/\*\*\*/,$cdsseq_cdscount);
				
				print "this cds sequence count is =", $count, "\n";
				print "this protein's cds seq is:\n",$cdssequence, "\n";

				creatMetCDSseqsTable($protein_id, $cdssequence );

				@prot_names = split (/\s/,$line);
				$protein_name = @prot_names[0];
				$protein_name =~ s/\>//g;

				$pep_signal = 0;
				$pepseq='';
				$protein_count+= 1;				
				print "protein count after increment = ", $protein_count, "\n";		       	
				next;
			    }		    
			    elsif(($line !~ /^\>/) && ($line =~ /^[A-Z]+$|^[A-Z]+\*$/)){
				#print "line =* ", $line, "\n";
				$pepseq .= $line;
				$pepseq =~ s/''|\W|\n|\s//g;
				$pep_signal = 1;
				
			    }

			} # 11 end of while loop

			print "Outside the loop, values for last protein are:\n",
			"protein name = ", $protein_name, "\n",
			"contig name = ", $contig_name, "\n",
			"protein sequence:\n",
			$pepseq, "\n";
			#$protein_count+= 1;

			my $length = @contigcdsseqs;
			print "total cds count of the contig = ", $length, "\n";
		       
			print "Total protein count of the contigs = ", $protein_count, "\n";

     			my $protein_id = creatMetProteinTable($pepseq, $contig_name, $protein_name);
			$cds = @contigcdsseqs[$protein_count-1];
			
			($cdssize, $cdsseq_cdscount)  = split(/:::/, $cds);
			($cdssequence, $count) = split(/\*\*\*/,$cdsseq_cdscount);
			
			print "this cds seq count is = ", $count, "\n";
			print "this protein's cds seq is:\n",$cdssequence, "\n";
			creatMetCDSseqsTable($protein_id,$cdssequence);

			
			print "insert value for total proteins column in protein table: \n", $protein_count,"\n",$contig_name, "\n"; 
			$query_handle = $connect_to_host->prepare('UPDATE contigs  SET proteins = ? WHERE name = ?');	
			$query_handle->execute($protein_count,$contig_name ) or die "could not execute statement: ".$query_handle;		 
			$protein_count = 0;
			@contigcdsseqs=();
			$pep_begin_signal = 0;
			$pep_signal = 0;
			$pepseq='';
			
		    }# 10

		}# 9
	    }# 8
	}# 7
	

	
	print "\n\n";
    }# 2
    
}# 1

close (READFILE);
	
	
    }# END OF METAGENOMICS DB
    
    $query_handle->finish;
    
    
}
    
    
    exit;



##################################################################################################


sub contigSeq{# 1
    my($file) = @_;
    my $contig_start_switch = 1;
    my $cd_range;
    my $tln_start = 0;
    my $tln_end = 0;
    my $protein_sequence;
    my $cd_count = 0;
    my $protein_count = 0;
    my $full_protein;
    my $this_protein_cdsrange;
    my $all_protein_set;
    my $contig_total_protein = 0;
    my @cds_range_array=();
    my $contig_name;
    my $contig_species;
    my $contig_taxonomy;
    my $contig_taxon_id;
    my $gene_signal = 0;
    my $org_signal;
    my $gene_name_duplicate;
    my $gene_name;
    my %hash_all_genes_cds_proteins_a_contig = ();
    my $cd_signal = 0;
    my $cd_product_signal = 0;
    my $protein_annotation;
    my $db_xref_hash_key;
    my $db_xref ;
    my $prot_flag;
    my $translation;
    my $cd_progress_count = 0;
    my $origin_switch = 0;
    my @contiglines =();
    my $contig_gc = 0;
    my %contig_all_features_hash=();
    my $contig_field_form= 'chromosome';
    my $xdb_signal = 0;
    my $comp_join_signal = 0;

    open (FILEHANDLE, "<$file") or die "can not open the file in sub : ", $file;
    while ( my $line = <FILEHANDLE>){ # 2

	next if ($line =~ /^\s+$/);
#	next if ($line =~ m/^\s+$|(misc_feature)|codon|note="|\s*protein_id|(="GeneID)|^\s+mRNA / );
	next if ($line =~ m/^\s+$|(misc_feature)|note="|\s*protein_id|(="GeneID)/ );
	    
	if (($line =~ /^VERSION/) && ($contig_start_switch == 1 )){ 
	    
	    ($contig_name )=  ($line =~ /^VERSION\s+(\S+)\s/);
	    print "New contig is being read ..............","\n";
	    print "contig_name = ",$contig_name, "\n"; 
	    
	    next;
	    
	}elsif (($line =~ /^SOURCE/)&& ($contig_start_switch == 1 )){
		
	    ($contig_species) =  ($line =~ /^SOURCE\s+(\S+.*$)/);
	    print "contig species = ", $contig_species, "\n";
		
	}elsif (($line =~ /ORGANISM/)&& ($contig_start_switch == 1 )){
	    $org_signal = 1;
	    
	    next;
	    
	}elsif(($org_signal ==1) && ($line !~ /REFERENCE/)){
	    $line =~ s/^\s+|\n$//;
	    
#	    $contig_taxonomy = $line;
	    push(@taxonomy, $line);
	    
	}elsif(($org_signal == 1) && ($line =~ /REFERENCE/)) {
	   
	    $contig_taxonomy = join('',@taxonomy);
 
	    $contig_taxonomy =~ s/\n//g;
	    $contig_taxonomy =~ s/\.$//;

#	    ($contig_taxonomy  = $contig_taxonomy) =~ s/\.$/$1/;
#	    $contig_taxonomy =~ s/\n//g;
	    print "taxonomy = ", $contig_taxonomy, "\n";	         
	    $org_signal = 0; 
	    @taxonomy =();
	    next;
	    
	} elsif($line =~ m/^\s*source/  ){
	    
	    $taxon_id_signal = 1;
	    
	}elsif (($line =~ m/^\s*\/db_xref\=\"(taxon)\:/i)&& ($contig_start_switch == 1 )){
	    
	    ($contig_taxon_id) = ($line =~ /taxon:(\d+)\"/);
	    print "taxon_id = ", $contig_taxon_id, "\n";
	    
	    next;
	}elsif (($line =~ /^\s+(gene)\s+/) && ($gene_signal == 0)){
	    
	    $gene_signal = 1; 
	    next;
# /locus_tag="	    
	}elsif ((($line =~ /\s+(\/gene=\")\S+.*\"$|\s+(\/locus_tag=")\S+.*\"$/) && ($gene_signal ==1  ))&& ($contig_start_switch == 1)){
	    
#	    ($gene_name) = ($line =~ /^\s+\/gene=\"(\S+.*)\"$/);
	    ($gene_name) = ($line =~ /^\s+\/\S.*=\"(\S+.*)\"$/);

	    next if($gene_name eq $gene_name_duplicate);
	    $gene_name_duplicate = $gene_name;
	    
	    print "\n#################\ngene name = ", $gene_name, "\n";
	    # $gene_name = $gene_name."@@@";
	    
	    $all_protein_set .= $gene_name."@@@";
	    $gene_signal = 0;
		
	    $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
		
#	}elsif (($gene_signal == 1) && ( $line !~ /\s+(\/gene=\")\S+.*\"$|\s+(\/locus_tag=")\S+.*\"$/ )) {
	}elsif (($gene_signal == 1) && ( $line !~ /\s+(\/gene=\")\S+.*\"$/ )) {
	    $gene_name = 'MISSING';
	    print "\ngene name = ", $gene_name, "\n";
	    #$gene_name = $gene_name."@@@";
	    
	    $all_protein_set .= $gene_name."@@@";
	    $gene_signal = 0;
	    
	    $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
	    
	}#elsif (($line =~ /(^\s*CDS\s*)/)&& ($contig_start_switch == 1 )){
	   # $line =~ s/\<|\>//g;
	   # $cd_signal = 1;
	   # $cd_product_signal = 1;

#DETERMINING IF THE CDS IS COMPLEMENT OR NOT AND GETTING THE CORRESPONDING SEQUENCE FROM CONTIG SEQUENCE	    
#	elsif (($line =~ /^\s+CDS\s+(complement)/ ) && ($contig_start_switch == 1 )){
	elsif (($line =~ /^\s+CDS\s+(complement)/ ) && ($line !~ /(join)/ ) ){
	    $cd_progress_count += 1;
	    $line =~ s/\<|\>//g;
	    $cd_signal = 1;
	    $cd_product_signal = 1;
	    $xdb_signal = 1;
		    
	    print "---- CD # $cd_progress_count.... is being read......\n";
	    
	    $line =~ s/\s+//g; 
	    ($cd_range) =	($line =~ /(\d+\.\.\d+)/);
	    
	    $cd_range = $cd_range."***comp";
	    print "cd_range from comp = ", $cd_range, "\n";
	    $this_protein_cdsrange = $cd_range;
	    
#		    push (@cds_range_array, $cd_range);
#		    $total_db_cds_count  += 1;
#		    print "size of cds range array for comp cds = ", my $length = @cds_range_array, "\n";
		   # $cd_range = $cd_range;
	    $all_protein_set .= $cd_range;
	    
	    $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
	     print  "cd_product_signal in cds comp line 7 = ", $cd_product_signal, "\n";
	    next;
	    
#	}

	}
	elsif ($line =~ /^\s+CDS\s+(complement\(join\(\d+)/  ){
	    $line =~ s/\<|\>//g;
	    $comp_join_signal = 1;
	    $cd_product_signal = 1;
	    $xdb_signal = 1;
	    push(@compjoinlines, $line);
	    print  "cd_product_signal in cds comp join first line 8 = ", $cd_product_signal, "\n";
	    next;
	    #}
	    
	}elsif(($comp_join_signal ==1) && ($line !~ /(^\s+\/\S.*$)/)){
	    $line =~ s/\<|\>//g;
	     print  "cd_product_signal in cds comp join continuation line 9 = ", $cd_product_signal, "\n";
	     push(@compjoinlines, $line);
	    next;
	    
	}elsif(($comp_join_signal ==1) && ($line =~ /(^\s+\/\S.*$)/)){
	    $comp_join_signal = 0;
	    $cd_signal = 1;
	    $cd_product_signal = 1;
	    $cd_progress_count += 1;
	    $comp_join_cdsranges = join ('',@compjoinlines);
	    
	    print "---- CD # $cd_progress_count.... is being read......\n";	

	    ($cd_range = $comp_join_cdsranges) =~ s/complement|CDS|join|\s+|\(|\)|\n|\<|\>/$1/g;
	    $cd_range = $cd_range."***comp+++join";    
	    print "the comp join cds ranges are :\n",
	    $cd_range, "\n";
	    
	    @compjoinlines=();
	    
	
	    $this_protein_cdsrange = $cd_range;
	    $all_protein_set .= $cd_range;
	    $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
	   print  "cd_product_signal in cds final comp join line 10 = ", $cd_product_signal, "\n";
	    
	    next;
	}


#	elsif (($line =~ /^\s+(CDS)\s+(join)\(\d+/ ) && ($contig_start_switch == 1 )){
	elsif ($line =~ /^\s+(CDS)\s+(join)\(\d+/ ){
	    $join_signal = 1;
	    $cd_product_signal = 1;
	    $xdb_signal = 1;
	    push(@joinlines, $line);
	    print  "cd_product_signal in cds join first line 6 = ", $cd_product_signal, "\n";
	    next;
	    #}
	    
	}elsif(($join_signal ==1) && ($line !~ /(^\s+\/\S.*$)/)){
	    print  "cd_product_signal in cds join continuation line 5 = ", $cd_product_signal, "\n";
	    push(@joinlines, $line);
	    next;
	    
	}elsif(($join_signal ==1) && ($line =~ /(^\s+\/\S.*$)/)){
	    $join_signal = 0;
	    $cd_signal = 1;
	    $cd_product_signal = 1;
	    $cd_progress_count += 1;
	    $join_cdsranges = join ('',@joinlines);
	    
	    print "---- CD # $cd_progress_count.... is being read......\n";	

	    ($cd_range = $join_cdsranges) =~ s/CDS|join|\s+|\(|\)|\n|\<|\>/$1/g;
	    $cd_range = $cd_range."***forward+++join";    
	    print "the join cds ranges are :\n",
	    $cd_range, "\n";
	    
	    @joinlines=();
	    
	
	    $this_protein_cdsrange = $cd_range;
	    $all_protein_set .= $cd_range;
	    $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
	    print  "cd_product_signal in cds final join line 4 = ", $cd_product_signal, "\n";
	    
	    next;
	}

#	elsif (($line =~ /^\s+(CDS)\s+(\d+\.\.\d+$)/ ) && ($line !~ /^\s+(CDS)\s+(join)\(\d+|^\s+(CDS)\s+(complement)\(\d+\.\.\d+\)$/)){
	    elsif (($line =~ /^\s+(CDS)\s+(\d+\.\.\d+$)|^\s+(CDS)\s+(\<\d+\.\.\d+$)|^\s+(CDS)\s+(\d+\.\.\d+\>$)/ ) && ($line !~ /^\s+(CDS)\s+(join)\(\d+|^\s+(CDS)\s+(complement)\(\d+\.\.\d+\)$/)){
		
	    $line =~ s/\<|\>//g;
	    $cd_signal = 1;
	    $cd_product_signal = 1;
	    $cd_progress_count += 1;
	    $xdb_signal = 1;
	   
	    print "---- CD # $cd_progress_count.... is being read......\n";	
	    
#EXTRACTING FORWARD STRAND CD INFORMATION FOR THE GENOME SEQUENCE
	    $line =~ s/\s+//g; 
	    
	    $line =~ /(\d+\.\.\d+)/;
	    $cd_range = $1;		
	    
	    $cd_range = $cd_range."+++positive***forward";
	    
	    print "cd_range from forward strand = ", $cd_range, "\n";
	    
	    $this_protein_cdsrange = $cd_range;
	    
#		push (@cds_range_array, $cd_range);
#		$total_db_cds_count  += 1;
#		print "length of cds range array = ", my $length = @cds_range_array, "\n";
	    #	$cd_range = $cd_range.":::";
		
	    $all_protein_set .= $cd_range;
	    $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
	    print "cd_product_signal in forward strand = ", $cd_product_signal, "\n";
	    next;
	}

#-----
	
	elsif (($line =~ /^(\s+\/product\=\"\S+.*\"$)/ ) && ($cd_product_signal == 1 )){
	    print "cd_product_signal in product line 1 = ", $cd_product_signal, "\n";
	    $cd_product_signal = 0;
	   # $cd_product_nextline_signal = 1;
	    $line =~ s/^\s+//;
	    $xdb_signal = 1;
	   # $pr_end_signal = 1;

	    ($protein_annotation ) = $line =~ /product\=\"(\S+.*)\"$/;
	  #  print "line is = ", $line, "\n";
#	    push(@cdproductlines, $line);
	    print "protein annotation from single product line = ",$protein_annotation, "\n";
	    next;

	}
#***#***

	elsif (($line =~ /^(\s+\/product\=\"\S+.*$)/ ) && ($cd_product_signal == 1 )){
	    print  "cd_product_signal in product line 2 = ", $cd_product_signal, "\n";
	    $cd_product_signal = 0;
	    $cd_product_nextline_signal = 1;
	    $line =~ s/^\s+//;
	    $xdb_signal = 1;
	    
	    ($line) = $line =~ /product\=\"(\S+.*$)/; 
	    $protein_annotation = $line;
	    push(@cdproductlines, $line);
	    print "collecting protein annotation from first line.... ", $line, "\n";
	    next;
	    
	}elsif (($cd_product_nextline_signal == 1 ) && ($line =~ /^(\s+\S+.*\"$)/  )){
	    print "cd_product_signal in product line 3 = ", $cd_product_signal, "\n";
	    $cd_product_nextline_signal = 0;
	    $cd_product_signal = 0;
	    $line =~ s/^(\s+)//;
	    $xdb_signal = 1;
	 
	  #  $protein_annotation = join ('',@cdproductlines);
	    $protein_annotation .= " ".$line;
	    ($protein_annotation = $protein_annotation) =~ s/"$/$1/g;
	    print "protein annotation from multilines = ",$protein_annotation, "\n";
	    @cdproductlines=();
	}
#***#***

#CAPTURING db_xref == NAME FEILD IN PROTEIN TABLE
###	elsif (($line =~ /(db_xref="GI:)/ )&& ($contig_start_switch == 1 )){
	elsif (($line =~ /(db_xref="GI:)/ )&& ($xdb_signal == 1 )){
	    $xdb_signal = 0;
	    $line =~ /db_xref\="GI:(\S+.*)"$/;
	    $db_xref = $1;
	    $db_xref=~ s/"//;
	    
	    print "protein name = ", $db_xref, "\n";
	    
	    $db_xref_hash_key = $db_xref;
	    
	    $db_xref = $db_xref."^^^";
	    $prot_flag = ":::".$db_xref.$protein_annotation;
	    
	    $all_protein_set .= ":::".$db_xref;

	    $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
	    
	    next;
	    
	}



# CAPTURING TRANSLATION == PROTEIN SEQ IN PROTEIN TABLE
# PROTEIN LINE BEGINS
#	elsif ($line =~ m/^\s+\/translation="([A-Z]+)$/){
	elsif (($line =~ m/^\s+\/(translation\=\"[A-Z]+$)/) && ($tln_start == 0 )){
#	    $db_protein_count += 1;
#	    $protein_count += 1;
	    $tln_start = 1;
	    $tln_end = 0;
	    
	    $line =~ /^\s+\/translation\=\"([A-Z]+)$/;	    
	    $translation = $1;
	    
	    $full_protein = $translation;
	    $protein_sequence = $translation;
	    
#	    print "protein sequence:\n", $protein_sequence, "\n";

	    next;
	    
#    DETECTING PROTEIN SEQ ENDED IN ONE LINE	    
#	}elsif ($line =~ m/^\s+\/(translation="[A-Z]+")$/){
#       }elsif ($line =~ m/^\s+\/(translation="[A-Z]+")$/){ # 07/10/08
	}
	
	elsif (($line =~ m/^\s+\/(translation\=\"[A-Z]+\")$/) && ($tln_start == 0 )){
	    
	    $line =~ /translation="([A-Z]+)\"$/;
	    $translation = $1;
	    
	    $full_protein = $translation;	    
	    
	    if(($cd_range !~ /(xxxxxx)/) && ($full_protein !~ /(xxxxxx)/) ) {
		$cd_count += 1;
		$total_db_cds_count  += 1;
		push (@cds_range_array, $cd_range);
		$db_protein_count += 1;
		$protein_count += 1;
		$tln_start = 0;
		$tln_end = 1;
		$contig_total_protein = $protein_count;

		$protein_sequence = $gene_name."@@@".$protein_count."===".$translation.$prot_flag."---".$cd_range;

#		$hash_get_protein{$this_protein_cdsrange} = $gene_name."@@@".$protein_count."===".$translation."###".$prot_flag;		
		$hash_get_protein{$this_protein_cdsrange} = $protein_sequence;

#		$protein_sequence = $protein_count."===".$translation."###";
		print "size of cds range array  = ", my $length = @cds_range_array, "\n";
		print "protein sequence:\n", $protein_sequence, "\n";
		$all_protein_set .= $protein_sequence;
		$hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
		
	    }
	    next;
	    
#SELLECTING PROTEIN SEQ MULTI LINES	    
#	}elsif ((($tln_start == 1) && ($tln_end == 0)) && ($line !~ /[A-Z]+"$/)) {
#	}elsif ((($tln_start == 1) && ($tln_end == 0)) && (($line =~ m/^\s+([A-Z]+)$/ )&& ($line !~ /[A-Z]+"$/))) { # 07/07/08
	    
	    
# DETECTING TRANSLATION CONTINUATION LINE
	    
	}
	
	elsif ((($tln_start == 1) && ($tln_end == 0)) && (($line =~ m/^\s+([A-Z]+$)/ )&& ($line !~ m/^\s+([A-Z]+\"$)/))) { # 07/10/08    
	    next if ($line =~ /codon|start|origin|gene|complement|\//);
	    
	    $line =~ /([A-Z]+$)/;	    
	    $translation = $1;
	    
#	    print $translation. "\n";
	    $full_protein .= $translation;
	    $protein_sequence .= $translation;
	    
#print "protein sequence:\n", $protein_sequence, "\n";
	    next;
	    
#DETECTING END OF PROTEIN SEQ LINE	    
#	}elsif (($line =~ /([A-Z]+)"$/) && ( $tln_end == 0)){ 07/10/08
	    
	}
	
	elsif (($line =~ m/(^\s+[A-Z]+\"$)/) && ( $tln_end == 0)){
	    
	    next if ($line =~ /codon|start|gene|complement|translation|db_xref|product|protein|locus|inference|organism|type|strain|\d\/|details/);
	    
#	    $line =~ /([A-Z]+)"$/; #07/10/08
	    $line =~ /^\s+([A-Z]+)\"$/;
	    
	    $translation = $1;
	    $protein_sequence .= $translation;
	    $full_protein .= $translation;

#	    $tln_end = 1; $tln_start = 0;

	    if(($cd_range !~ /(xxxxxx)/) && ($full_protein !~ /(xxxxxx)/) ) {
		$cd_count += 1;
		$total_db_cds_count  += 1;
		push (@cds_range_array, $cd_range);
		$db_protein_count += 1;
		$protein_count += 1;
		$tln_start = 0;
		$tln_end = 1;
		$contig_total_protein = $protein_count;
       

		$protein_sequence = $gene_name."@@@".$protein_count."===".$protein_sequence.$prot_flag."---".$cd_range;

#		$hash_get_protein{$this_protein_cdsrange} = $gene_name."@@@".$protein_count."===".$translation.$prot_flag;		
		$hash_get_protein{$this_protein_cdsrange} = $protein_sequence;

#		$protein_sequence .= $translation."###"; 
#		$protein_sequence = $protein_count."===".$protein_sequence;
		print "size of cds range array  = ", my $length = @cds_range_array, "\n";
		print "protein sequence:\n", $protein_sequence, "\n";
		print "contig start switch = ", $contig_start_switch, "\n";
		$all_protein_set .= $protein_sequence;
		$hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
	       
		next;
	    }



#	    $hash_get_protein{$this_protein_cdsrange} = $gene_name."@@@".$protein_sequence.$translation.$prot_flag;
	    
#	    $protein_sequence .= $translation."###";    
	    
#	    print "protein sequence:\n", $protein_sequence, "\n";

#	    $all_protein_set .= $protein_sequence;
#	    $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
#	    $contig_total_protein = $protein_count;

	    
	    
	}#elsif (($tln_end == 1) && ($tln_start == 0) ){
	    
	  
	  #  $tln_start = 0;
	  #  $tln_end = 0;
	    
	    # $hash_all_genes_cds_proteins_a_contig{$contig_name} = $all_protein_set;
	    
#	}
	
#***	    
	elsif (($line =~ /^ORIGIN/ )&& ($contig_start_switch == 1)) {
	    
	    print "\nI am at the origin line.... $line", "\n";
	    if (($cd_count == 0) || ($protein_count == 0)){
	 
	    $cd_range = 'xxxxxx';
	    $tln_start = 0;
	    $tln_end = 0;
	    $protein_sequence='xxxxxx';
	    $full_protein ='xxxxxx';
	    $this_protein_cdsrange='';
	    $all_protein_set='';
	    $contig_name='xxxxxx';
	    $contig_species='xxxxxx';
	    $contig_taxonomy='xxxxxx';
	    $contig_taxon_id='xxxxxx';
	    $gene_signal = 0;
	    $org_signal;
	    $gene_name_duplicate;
	    $gene_name='xxxxxx';
	    $cd_signal = 0;
	    $cd_product_signal = 0;
	    $protein_annotation='xxxxxx';
	    $db_xref_hash_key='';
	    $db_xref='xxxxxx';
	    $prot_flag='';
	    $translation='';
	    $contig_gc=0;
	    @contiglines = ();
	    $contig_table_insert_values='xxxxxx';
	    %contig_all_features_hash=();

	    next;
	    }
	    #   $contig_start_switch = 0;
	    $cd_range = 'xxxxxx';
	    $tln_start = 0;
	    $tln_end = 0;
	    $this_protein_cdsrange='';
	    $all_protein_set='';
	    $gene_signal = 0;
	    $org_signal;
	    $gene_name_duplicate;
	    $cd_signal = 0;
	    $cd_product_signal = 0;
	    $db_xref_hash_key='';
	    $prot_flag='';
	    $translation='';
	    $origin_switch = 1;
	    $contig_line ='';   
	    $contig_seq = 'xxxxxx';
	    @contiglines = ();
	    $origin_switch = 1;
	    $contig_gc=0;
	    $contig_table_insert_values='xxxxxx';
	    %contig_all_features_hash=();

		
	}elsif( ($line !~ /^\/\/\n/  )&& ($origin_switch == 1)){
	   	   
	   $line  =~ s/\s//g;
	   $line =~ /\d+([a-z]+)$/;
	   $contig_line = $1;
	   $contig_line =~ s/\n//g;	   
	   push (@contiglines, $contig_line);

       
	   $contig_start_switch = 0;
	}
	elsif(($line =~ /^\/\/\n/  ) && ($contig_start_switch == 0 ) ){
	    $origin_switch = 0;
	    $genome_count += 1;
#process contig's all info

	    $contig_seq = join('', @contiglines);
	    $contig_seq =~ s/\W|\s|\n|\d//g;
	    $contig_size = length $contig_seq;;
	    print "Lenth of the contig is = ", $contig_size, "\n";
#****	    print "Contig Sequence:", "\n", $contig_seq, "\n";
 
#CALCULATING GC CONTENT OF CONTIG
	    $g = ($contig_seq =~ tr/Gg//);
	    $c = ($contig_seq =~ tr/Cc//);
	    $contig_gc = $g + $c;
	    $contig_gc = ($contig_gc/length $contig_seq)*100;

	    $contig_table_insert_values = $contig_size.":::".$contig_species.":::".$contig_taxon_id.":::".$contig_taxonomy.":::".$contig_field_form.":::".$contig_total_protein.":::".$contig_gc;
	    
	    $contig_all_features_hash{$contig_name } = $contig_table_insert_values;

	    print "I am at the contig's end //line...processing its information\n";
	    print "This Contig's total proteins = ", $contig_total_protein, "\n";
	    print "This contig's total cds seqs counted in protein section= ", $cd_count, "\n";
	    print "This contig's total cds counted while reading each cds seq = ", $cd_progress_count, "\n";
#****	    print "contig insert values :\n", $contig_table_insert_values, "\n";
	    print "Total proteins in database counted = ", $db_protein_count, "\n";
	    print "total database CDS seqs counted = ", $total_db_cds_count, "\n";
	    print "size of contig's cds range array = ", my $length = @cds_range_array, "\n";
	    print "CONTIG NUMBER = $genome_count", " HAS BEEN PROCECESSED......\n";

# CREATING TABLE
	    $table_status = creatContigTable(\%contig_all_features_hash, $contig_name, $contig_seq);
	    print $table_status, "\n";	    

# CREATING MOLECULE TABLE
	    $molecule_table_status = createMoleculesTable ($contig_seq, $contig_name );
	    print $molecule_table_status;
	    
# PROCESSING CDS RANGES
	      processCDS(\@cds_range_array, \%hash_get_protein, $contig_seq,$contig_name,\%hash_all_genes_cds_proteins_a_contig  );





# REGENEARTING VARIABLES------------------------------------------------------------------	    
	    $contig_start_switch = 1;
	    $protein_sequence = 'xxxxxx';
	    $cd_range = 'xxxxxx';
	    $tln_start = 0;
	    $tln_end = 0;
	    $cd_count = 0;
	    $protein_count = 0;
	    $full_protein = 'xxxxxx';
	    $this_protein_cdsrange='';
	    $all_protein_set='';
	    $contig_total_protein = 0;
	    @cds_range_array=();
	    $contig_name='xxxxxx';
	    $contig_species='xxxxxx';
	    $contig_taxonomy='xxxxxx';
	    $contig_taxon_id='xxxxxx';
	    $db_xref ='xxxxxx';
	    $gene_name = 'xxxxxx';
	    $hash_all_genes_cds_proteins_a_contig = ();
	    $cd_signal = 0;
	    $cd_product_signal = 0;
	    $protein_annotation = 'xxxxxx';
	    $db_xref_hash_key='';
	    $prot_flag ='';
	    $translation='';
	    $cd_progress_count = 0;
	    $contig_line='';   
	    $contig_seq = 'xxxxxx';
	    @contiglines=();
	    $contig_gc=0;
	    $contig_table_insert_values='xxxxxx';
	    %contig_all_features_hash=();
	    $xdb_signal = 0;

	    print "The size of cds range array after contig processing is = ", my $length = @cds_range_array, "\n\n";
		

	}
	
	
    }# 2 END Of WHILE LOOP CALLING FILEHANDLE    
    
    
    return ($genome_count, $total_db_cds_count, $db_protein_count );

}# 1 End of sub contigSeq($name){



sub creatContigTable{

    my ($contig_features_ref, $contig_name_ref, $contig_seq_ref) = @_;
    
    my  %contig_table_fields_hash = %$contig_features_ref;
    my  $contig_name = $contig_name_ref;
    my $contig_seq = $contig_seq_ref;
    
    print "contig name in the sub create table = ", $contig_name, "\n";

    my $insert_values =  $contig_table_fields_hash{$contig_name};

    print "insert values in create table sub are:\n",$insert_values, "\n";

# $contig_table_insert_values = $contig_size.":::".$contig_species.":::".$contig_taxon_id.":::".$contig_taxonomy.":::".$contig_field_form.":::".$contig_total_protein.":::".$contig_gc;

 my   ($contig_size, $contig_species, $contig_taxon_id, $contig_taxonomy, $contig_field_form, $contig_total_protein, $contig_gc) = split (/:::/, $insert_values);
    

# INSERTING....    
    $query_handle = $connect_to_host->prepare('INSERT into contigs (name, size, species, taxon_id, taxonomy, form, proteins, gc)VALUES (?,?,?,?,?,?,?,?)');
   $query_handle->execute($contig_name, $contig_size, $contig_species, $contig_taxon_id, $contig_taxonomy, $contig_field_form, $contig_total_protein, $contig_gc ) or die "could not execute statement: ".$query_handle;
    
    print "insert IN CONTIGS successful", "\n";




#    $query_handle->finish;
#FINISHED INSERTING    

#    createMoleculesTable ($contig_seq, $contig_name );

    return "Contig info in contig table inserted successfully\n";
 
}# END OF CONTIG TABLE SUB


sub createMoleculesTable{
    my ($contig_seq_ref,$contig_name_ref) = @_;
    my $contig_seq = $contig_seq_ref;
    my $contig_name = $contig_name_ref;
    
    print "contig name in create molecules table = ", $contig_name, "\n";

# my $sth = $dbh->prepare('SELECT * FROM people WHERE lastname = ?')
#$query = "SELECT id, name FROM contigs";
    
    $query_handle = $connect_to_host->prepare('SELECT id, name FROM contigs WHERE name = ?');
    $query_handle->execute($contig_name) or die "could not execute statement in sub createMoleculesTable: ".$query_handle;
    @data = $query_handle->fetchrow_array();
	
    print "The data extracted from contig's table are:\n";
    print @data[0]."=".@data[1],"\n"; # @data[1]= name of contig, @data[0]= contig_id from contig table
    
    $query_handle = $connect_to_host->prepare('INSERT INTO molecules  (contig_id, seq)VALUES (?,? )');
    $query_handle->execute(@data[0],$contig_seq) or die "could not execute statement in sub createMoleculesTable: ".$query_handle;
    
    print " I have successfully inserted contig's sequence in molecules table in the sub createMoleculesTable \n";
    
    return "I have successfully inserted contig's sequence in molecules table in the sub createMoleculesTable\n"

}# END OF sub createMoleculesTable{



sub processCDS{
   
    my ($array_ref, $hash_ref, $contigseq_ref, $contig_name_ref, $hash_all_genes_cds_proteins_a_contig_ref) = @_;    
    my @cds_range_array_sub = @$array_ref;
    my %get_protein_hash_sub = %$hash_ref;
    my $contig_seq = $contigseq_ref;
    my $contig_name = $contig_name_ref;
    my %hash_all_genes_cds_proteins_a_contig = %$hash_all_genes_cds_proteins_a_contig_ref;
    my ($contig_id, $protein_id, $protannotn)='';


#   my $prot_info =  $hash_all_genes_cds_proteins_a_contig{$contig_name};
    
#    print "the protien info present in hash in sub processCDS :\n$prot_info", "\n";

    print "\nThe length of contig in sub processCDS = ", length $contig_seq, "\n";

    my $cd_count = 0;
    print "the length of the cds range array in sub processCDS = ", $length = @cds_range_array_sub, "\n";
 

    foreach my $cds_range (@cds_range_array_sub){
	$cd_count += 1;
	my $prot_seq_string = $get_protein_hash_sub{$cds_range};

	my ($gene_nm_prot_num_prot_seq, $prot_name_prot_annot_cdsrange_string)  = split(/:::/, $prot_seq_string);

	my ($gene_nm, $prot_num_prot_seq) = split (/\@\@\@/, $gene_nm_prot_num_prot_seq ); 

	print "\nfor the cd range = $cds_range\nprotein sequence extracted is:\n$prot_num_prot_seq\n" ;
	
	print "\n######## FROM CDRANGE ARRAY, CD number = ", $cd_count, " IS BEING PROCESSED NOW .....\n";
	
	if (($cds_range =~ m/\d+\.\.\d+\*\*\*(comp)$/ ) && ($cds_range !~ /(comp\+\+\+join)$/ )) {
	    
	    $cds_range =~ /(\d+\.\.\d+)\*\*\*comp$/;
	    my ($cd_range_extracted) = $1;
	    
#getting reverse comp cds sequence
	    my ($start_cds_stop_string) = rev_compSeq($cd_range_extracted, $contig_seq);   
	    
	    my ($start,$cds_rev_comp,$stop) = split (/:::/,$start_cds_stop_string);
	    my ($cds_seq) = $cds_rev_comp;
	    
	    
#PUTTING TOGETHER CDS RANGE to USE IN CDS TABLE		       
	    my  $this_cds_range = $start."..".$stop;
	    
	    print "****BEGINS******printing from sub processCDS, this is COMP cds, \nits start position = ", $start, "\n", 
	    "cds sequence to translate: \n", $cds_seq, "\n", "the stop position = ", $stop, "\n";
	    print "So, the cds range is : ", $this_cds_range, "\n\n";
	    
	    print "Length of cds sequence is = ", length $cds_seq, "\n";
	    
	    $cds_seq =~ s/\s|\n|\W//g;

	    print "the cds sequence has been cleaned and is being used for translation...\n";
	    
	    
#CURRENT TRANSLATION
	    $current_translated_protein = translateProtein ($cds_seq);
	    $current_translated_protein =~ s/\_$|\s|\W\n//;
	    
#	       my $protein_set = $hash_all_genes_cds_proteins_a_contig{$contig_name};
	    
	    my ($prot_number, $protein_sequence_extracted) = split (/===/,$prot_num_prot_seq);
	    
	    print "\n!!!!! I N ***** C O M P ***** CDS, Extracted protein:\n\n*", $protein_sequence_extracted, "\n-";
	    print $current_translated_protein, "\n\n", "The above was translated protein\n\n";
	    
	    if ($current_translated_protein eq $protein_sequence_extracted){
		
		print " ********* !!!  E X A C T - M A T C H !!! *******", "\n\n";
		
	    }else{
		
		print " ^^^^^^ N O T  S A M E ^^^^^", "\n\n";
		
	    }
	    

# CALLING FUNCTION createProteinsTable
#	    createProteinTable(\%hash_all_genes_cds_proteins_a_contig,$contig_name, $protein_sequence_extracted);
	   ($contig_id, $protein_id, $protannotn)= createProteinTable($contig_name, $protein_sequence_extracted, $prot_name_prot_annot_cdsrange_string);

	    print "From  processCDS sub, did I pass protein id value and cds seq to sub createCDSseqTable:\n",
	    "protein id = ", $protein_id, "\n",
	    "cds seq = ", $cds_seq, "\n\n";

# CALLING FUNCTION CDSseqTable
	   my $cdstable_status =  createCDSseqTable($protein_id, $cds_seq); 
	    print $cdstable_status, "\n";

# CALLING FUNCTION createGeneorderTable
	    if($gene_nm =~ /(MISSING)/){
		$gene_nm = $protannotn;
	    }

	    $start_cds_stop_string = $start_cds_stop_string."***-1";
	    createGeneorderTable($protein_id,$contig_id,$start_cds_stop_string,$gene_nm );

	    print "\n#### END ####CDS SEQUENCE AND PROTEIN SEQ INFO COLLECTION COMPLETED***** FROM \"comp\" CDS *******\n";
		   
	    next;
	    
	}#endof if statement used for checking cds range--- COMP
	
# CDS RANGE FOR comp+++join
	elsif (($cds_range =~ m/(comp\+\+\+join)$/) && ($cds_range !~ m/(comp)$|(forward\+\+\+join)$/)){
	    
	    $start_cds_stop_string = compJoinCds($cds_range, $contig_seq);
	    
	    ($start,$cds_seq,$stop) = split (/:::/,$start_cds_stop_string);
	    
	    my ($this_cds_range) = $start."..".$stop;
	    
#PUTTING TOGETHER CDS RANGE to USE IN CDS TABLE		       
	    print "****BEGINS******printing from sub processCDS, this is JOIN cds, \nits start position = ", $start, "\n", 
	    "cds sequence to translate: \n", $cds_seq, "\n", "the stop position = ", $stop, "\n";
	    print "So, the cds range is : ", $this_cds_range, "\n\n";
	    
	    print "Length of cds sequence is = ", length $cds_seq, "\n";
	    
	    $cds_seq =~ s/\s|\n|\W//g;
	    
	    print "the cds sequence has been cleaned and is being used for translation...\n";

	    
#CURRENT TRANSLATION
	    $current_translated_protein = translateProtein ($cds_seq);
	    $current_translated_protein =~ s/\_$//;
	    
#	       my $protein_set = $hash_all_genes_cds_proteins_a_contig{$contig_name};
	    
	    my ($prot_number, $protein_sequence_extracted) = split (/===/,$prot_num_prot_seq);
	    
	    print "\n!!!!! I N ***** J O I N  ***** CDS, Extracted protein:\n\n*", $protein_sequence_extracted, "\n-";
	    print $current_translated_protein, "\n\n", "The above was translated protein\n\n";
	    
	    if ($current_translated_protein eq $protein_sequence_extracted){
		
		print " ********* !!!  E X A C T - M A T C H !!! *******", "\n\n";
		
	    }else{
		
		print " ^^^^^^ N O T  S A M E ^^^^^", "\n\n";
		
	    }
	    
# CALLING FUNCTION createProteinsTable
#	    createProteinTable(\%hash_all_genes_cds_proteins_a_contig,$contig_name, $protein_sequence_extracted);
	    ($contig_id, $protein_id, $protannotn) = createProteinTable($contig_name, $protein_sequence_extracted, $prot_name_prot_annot_cdsrange_string);

# CALLING FUNCTION CDSseqTable
	   my $cdstable_status = createCDSseqTable($protein_id, $cds_seq);
	    print $cdstable_status, "\n";
	    if($gene_nm =~ /(MISSING)/){
		$gene_nm = $protannotn;
	    }

# CALLING FUNCTION createGeneorderTable
	    $start_cds_stop_string = $start_cds_stop_string."***-1";
	    createGeneorderTable($protein_id,$contig_id,$start_cds_stop_string, $gene_nm );

	    print "\n#### END ####CDS SEQUENCE AND PROTEIN SEQ INFO COLLECTION COMPLETED***** FROM \"join\" CDS *******\n";	       
	    
	    next;  
	    
#PROCESSING FORWARD CDS RANGES
	   }

# JOIN CDS RANGES
	elsif (($cds_range =~ m/(forward\+\+\+join)$/) && ($cds_range !~ m/(comp)$|(comp\+\+\+join)$/)){
	    
	    $start_cds_stop_string = joinCds($cds_range, $contig_seq);

	    ($start,$cds_seq,$stop) = split (/:::/,$start_cds_stop_string);
	    
	    my ($this_cds_range) = $start."..".$stop;
	    
#PUTTING TOGETHER CDS RANGE to USE IN CDS TABLE		       
	    print "****BEGINS******printing from sub processCDS, this is JOIN cds, \nits start position = ", $start, "\n", 
	    "cds sequence to translate: \n", $cds_seq, "\n", "the stop position = ", $stop, "\n";
	    print "So, the cds range is : ", $this_cds_range, "\n\n";
	    
	    print "Length of cds sequence is = ", length $cds_seq, "\n";
	    
	    $cds_seq =~ s/\s|\n|\W//g;
	    
	    print "the cds sequence has been cleaned and is being used for translation...\n";

	    
#CURRENT TRANSLATION
	    $current_translated_protein = translateProtein ($cds_seq);
	    $current_translated_protein =~ s/\_$//;
	    
#	       my $protein_set = $hash_all_genes_cds_proteins_a_contig{$contig_name};
	    
	    my ($prot_number, $protein_sequence_extracted) = split (/===/,$prot_num_prot_seq);
	    
	    print "\n!!!!! I N ***** J O I N  ***** CDS, Extracted protein:\n\n*", $protein_sequence_extracted, "\n-";
	    print $current_translated_protein, "\n\n", "The above was translated protein\n\n";
	    
	    if ($current_translated_protein eq $protein_sequence_extracted){
		
		print " ********* !!!  E X A C T - M A T C H !!! *******", "\n\n";
		
	    }else{
		    
		print " ^^^^^^ N O T  S A M E ^^^^^", "\n\n";
		
	    }
                   	    
# CALLING FUNCTION createProteinsTable
#	    createProteinTable(\%hash_all_genes_cds_proteins_a_contig,$contig_name, $protein_sequence_extracted);
	    ($contig_id, $protein_id, $protannotn) = createProteinTable($contig_name, $protein_sequence_extracted, $prot_name_prot_annot_cdsrange_string);

# CALLING FUNCTION createGeneorderTable
	    createCDSseqTable($protein_id, $cds_seq);

	    if($gene_nm =~ /(MISSING)/){
		$gene_nm = $protannotn;
	    }

# CALLING FUNCTION createGeneorderTable
	    $start_cds_stop_string = $start_cds_stop_string."***+1";
	    createGeneorderTable($protein_id,$contig_id,$start_cds_stop_string, $gene_nm );

	    print "\n#### END ####CDS SEQUENCE AND PROTEIN SEQ INFO COLLECTION COMPLETED***** FROM \"join\" CDS *******\n";	       
	    
	    next;  
	    
#PROCESSING FORWARD CDS RANGES
	   }

#	elsif($cd =~ m/(\d+\.\.\d+\*\*\*forward)$/  ){
	elsif (($cds_range =~ m/(positive)/) && ($cds_range !~ m/(comp)|(join)/)){
	    
	    $cds_range =~ /(\d+\.\.\d+)/;
	    my ( $cd_range_extracted) = $1;
	    
#geting cds seq
	    my $start_cds_stop_string = cdsSeq($cd_range_extracted, $contig_seq);
	    
	    ($start,$cds_seq,$stop) = split (/:::/,$start_cds_stop_string);
	    
#PUTTING TOGETHER CDS RANGE to USE IN CDS TABLE		       
	    my ($this_cds_range) = $start."..".$stop;
	    
#PUTTING TOGETHER CDS RANGE to USE IN CDS TABLE		       
	    print "****BEGINS******printing from sub processCDS, this is POSITIVE FORWARD cds, \nits start position = ", $start, "\n", 
	       "cds sequence to translate: \n", $cds_seq, "\n", "the stop position = ", $stop, "\n";
	    print "So, the cds range is : ", $this_cds_range, "\n\n";
	    
	    print "Length of cds sequence is = ", length $cds_seq, "\n";
	    
	    $cds_seq =~ s/\s|\n|\W//g;
	       
	    print "the cds sequence has been cleaned and is being used for translation...\n";
	    
#CURRENT TRANSLATION
	    $current_translated_protein = translateProtein ($cds_seq);
	    $current_translated_protein =~ s/\_$//;
	    
#	       my $protein_set = $hash_all_genes_cds_proteins_a_contig{$contig_name};
	    
	    my ($prot_number, $protein_sequence_extracted) = split (/===/,$prot_num_prot_seq);
	    
	    print "\n!!!!! I N ***** P O S I T I V E -  F O R W A R D  ***** CDS, Extracted protein:\n\n*", $protein_sequence_extracted, "\n-";
	       print $current_translated_protein, "\n\n", "The above was translated protein\n\n";
	    
	    if ($current_translated_protein eq $protein_sequence_extracted){
		
		print " ********* !!!  E X A C T - M A T C H !!! *******", "\n\n";
		
	    }else{
		
		print " ^^^^^^ N O T  S A M E ^^^^^", "\n\n";
		
	    }

 # CALLING FUNCTION createProteinsTable
#	    createProteinTable(\%hash_all_genes_cds_proteins_a_contig,$contig_name, $protein_sequence_extracted);
	   ($contig_id, $protein_id, $protannotn) = createProteinTable($contig_name, $protein_sequence_extracted, $prot_name_prot_annot_cdsrange_string);

# CALLING FUNCTION createGeneorderTable
	    my $cdstable_status = createCDSseqTable($protein_id, $cds_seq);
	    print $cdstable_status, "\n";

# CALLING FUNCTION createGeneorderTable
	    if($gene_nm =~ /(MISSING)/){
		$gene_nm = $protannotn;
	    }

	   $start_cds_stop_string = $start_cds_stop_string."***+1";
	    createGeneorderTable($protein_id,$contig_id,$start_cds_stop_string, $gene_nm );

	    print "\n#### END ####CDS SEQUENCE AND PROTEIN SEQ INFO COLLECTION COMPLETED***** FROM \"positive forward\" CDS *******\n";
	    
	    next;
	    
	}# end of if-else ($cd =~ m/\d+\.\.\d+\*\*\*(comp)$/ ){


    }# END OF cds_range array begin loop
    
}# end of sub processCDS




sub createProteinTable{

    my ($contig_name_ref, $protein_sequence_extracted_ref, $prot_annot_string_ref) = @_;

#    %hash_all_genes_cds_proteins_a_contig = %$hash_all_genes_cds_proteins_a_contig_ref;
    my  $prot_annot_string = $prot_annot_string_ref;
    my $contig_name = $contig_name_ref;
    my $protein_sequence = $protein_sequence_extracted_ref;

#    my $protein_info =  $hash_all_genes_cds_proteins_a_contig{$contig_name};
    
    my ($protein_name, $protein_annot_str ) = split (/\^\^\^/,$prot_annot_string);
    my ($protein_annotation, $cdrange) = split(/---/,$protein_annot_str);
     
    print "printing from sub createProteinTable:\n";
    print "Contig name is = ", $contig_name, "\n";
    print "protein sequence is : ", $protein_sequence, "\n";
    print "protein name = ",$protein_name, "\n";
    print "protein annotation = ", $protein_annotation, "\n";
 
# my $sth = $dbh->prepare('SELECT * FROM people WHERE lastname = ?')
#$query = "SELECT id, name FROM contigs";
    
    $query_handle = $connect_to_host->prepare('SELECT id, name FROM contigs WHERE name = ?');
    $query_handle->execute($contig_name) or die "could not execute statement in sub createProteinsTable: ".$query_handle;
    
    @data = $query_handle->fetchrow_array();
    
    print "The data extracted from contig's table in sub createProteinTable are:\n";
    print @data[0]."=".@data[1],"\n"; # @data[1]= name of contig, @data[0]= contig_id from contig table
    
    $query_handle = $connect_to_host->prepare('INSERT INTO proteins  (name, contig_id, annotation, seq)VALUES (?,?,?,? )');
    $query_handle->execute($protein_name,@data[0],$protein_annotation,$protein_sequence) or die "could not execute statement in sub createProteinTable: ".$query_handle;
    
    print " I have successfully inserted protein sequence info in the sub createProteinTable\n";

    
    $query_handle1 = $connect_to_host->prepare('SELECT id, name FROM proteins WHERE name = ?');
    $query_handle1->execute($protein_name) or die "could not execute statement in sub createProteinsTable: ".$query_handle1;
    
    @data_protein = $query_handle1->fetchrow_array();

    print "The protein data extracted from proteins table in sub createProteinTable are:\n";
    print "protein id = ", @data_protein[0]."="."protein name = ", @data_protein[1], "\n";

    $query_handle1->finish;
    return (@data[0], @data_protein[0], $protein_annotation);
    
}# END OF SUB createProteinTable



#CREATING CDS SEQ TABLE
sub createCDSseqTable{ #($protein_id, $cds_seq)

    my ($protein_id_ref, $cds_seq_ref) = @_;    
    my $protein_id = $protein_id_ref;
    my $cds_seq = $cds_seq_ref;

#INSERTING INTO CDSSEQS TABLE
        
    $query_handle = $connect_to_host->prepare('INSERT INTO  cdsseqs (protein_id, seq) VALUES (?,? )');
    $query_handle->execute($protein_id, $cds_seq ) or die "could not execute statementin sub createCDSseqTable : ".$query_handle;
      

    return "SUCCESSFULLY CREATED CDSSEQ TABLE\n";


} # endof SUB createCdseqTable


# createGeneorderTable($protein_id,$contig_id,$start_cds_stop_string, $gene_nm );
sub createGeneorderTable{

    my ($protein_id_ref,$contig_id_ref,$start_cds_stop_string_ref, $gene_nm_ref ) = @_;
    my $protein_id = $protein_id_ref;
    my $contig_id = $contig_id_ref;
    my $start_cds_stop_string_strand = $start_cds_stop_string_ref;
    my $gene_name = $gene_nm_ref;
    my $strand;
    
    ($start_cds_stop_string, $strand) = split (/\*\*\*/, $start_cds_stop_string_strand );
    
    my ($start,$cds_seq,$stop) = split (/:::/,$start_cds_stop_string);
    
    $query_handle = $connect_to_host->prepare('INSERT INTO  geneorders (protein_id, contig_id, start,stop,strand,gene_name )VALUES (?,?,?,?,?,? )');
    $query_handle->execute($protein_id,$contig_id,$start,$stop,$strand,$gene_name ) or die "could not execute statement in sub createCDSseqTable: ".$query_handle;
    
    print " I have successfully inserted geneorder info for the protein from the sub createGeneorderTable\n";

    $query_handle->finish;

}#END OF sub createGeneorderTable(




sub rev_compSeq{
	my ($cds, $cont_seq) = @_;
	my ($cds_sq);
#    print "cds  range = ", $cds, "\n";
#    print "contig in sub= ", $cont_seq, "\n";
	
	
	$cds =~ s/\s|\n|NULL//g;
	
	my ($cds_start, $cds_end) = split( /\.\./, $cds);
	
	$cds_sq  = substr ($cont_seq, $cds_start, $cds_end - $cds_start);
	
	$rev_cds_sq = reverse $cds_sq;
	
#    print "the rev cds sq in sub is = ", $rev_cds_sq, "\n";
	
	$rev_cds_sq =~ tr/ACGTacgt/TGCAtgca/;
	
#    print "rev comp seq in sub = ", $rev_cds_sq, "\n";
	
	return $cds_start.":::".$rev_cds_sq.":::".$cds_end;
	
	
}



sub translateProtein {
    
    my ($cds) = @_;
    my ($protein);
    my ($codon);
    
    $cds =~ s/\s|\n|NULL|\W//g;
    
    
#    print "cds length in tranlation sub = ", length $cds, "\n";
    
    for ($i = 0; $i < (length ($cds) - 2); $i+= 3){
#	print "i = ", $i, "\n";
	
	$codon = substr ($cds, $i, 3);
	$protein .=  codon2aa ($codon);
#	print "protein = ", $protein, "\n";
    }   
    
    return $protein;
}

sub codon2aa {
    my($codon) = @_;
    
    $codon = uc $codon;
    
    my(%genetic_code) = (
	
	'TCA' => 'S',    # Serine
	'TCC' => 'S',    # Serine
	'TCG' => 'S',    # Serine
	'TCT' => 'S',    # Serine
	'TTC' => 'F',    # Phenylalanine
	'TTT' => 'F',    # Phenylalanine
	'TTA' => 'L',    # Leucine
	'TTG' => 'L',    # Leucine
	'TAC' => 'Y',    # Tyrosine
	'TAT' => 'Y',    # Tyrosine
	'TAA' => '_',    # Stop
	'TAG' => '_',    # Stop
	'TGC' => 'C',    # Cysteine
	'TGT' => 'C',    # Cysteine
	'TGA' => '_',    # Stop
	'TGG' => 'W',    # Tryptophan
	'CTA' => 'L',    # Leucine
	'CTC' => 'L',    # Leucine
	'CTG' => 'L',    # Leucine
	'CTT' => 'L',    # Leucine
	'CCA' => 'P',    # Proline
	'CCC' => 'P',    # Proline
	'CCG' => 'P',    # Proline
	'CCT' => 'P',    # Proline
	'CAC' => 'H',    # Histidine
	'CAT' => 'H',    # Histidine
	'CAA' => 'Q',    # Glutamine
	'CAG' => 'Q',    # Glutamine
	'CGA' => 'R',    # Arginine
	'CGC' => 'R',    # Arginine
	'CGG' => 'R',    # Arginine
	'CGT' => 'R',    # Arginine
	'ATA' => 'I',    # Isoleucine
	'ATC' => 'I',    # Isoleucine
	'ATT' => 'I',    # Isoleucine
	'ATG' => 'M',    # Methionine
	'ACA' => 'T',    # Threonine
	'ACC' => 'T',    # Threonine
	'ACG' => 'T',    # Threonine
	'ACT' => 'T',    # Threonine
	'AAC' => 'N',    # Asparagine
	'AAT' => 'N',    # Asparagine
	'AAA' => 'K',    # Lysine
	'AAG' => 'K',    # Lysine
	'AGC' => 'S',    # Serine
	'AGT' => 'S',    # Serine
	'AGA' => 'R',    # Arginine
	'AGG' => 'R',    # Arginine
	'GTA' => 'V',    # Valine
	'GTC' => 'V',    # Valine
	'GTG' => 'V',    # Valine
	'GTT' => 'V',    # Valine
	'GCA' => 'A',    # Alanine
	'GCC' => 'A',    # Alanine
	'GCG' => 'A',    # Alanine
	'GCT' => 'A',    # Alanine
	'GAC' => 'D',    # Aspartic Acid
	'GAT' => 'D',    # Aspartic Acid
	'GAA' => 'E',    # Glutamic Acid
	'GAG' => 'E',    # Glutamic Acid
	'GGA' => 'G',    # Glycine
	'GGC' => 'G',    # Glycine
	'GGG' => 'G',    # Glycine
	'GGT' => 'G',    # Glycine
	);
    
    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
	print "##############^^^^^^^^^^^^^^^^^^^^*****************************\n";
#	print STDERR "Bad codon \"$codon\"!!\n";
            #exit;
    }
}

sub compJoinCds{

    my ($cd_ranges, $cont_seq) = @_;
    my $cds_sq;
    my $cd_start;
    my $cd_stop;
    
    ($cd_ranges_set  =  $cd_ranges) =~ s/\+|\*|join|comp/$1/g;
    
    my @array_cdranges = split (/\,/, $cd_ranges_set);

    print "printing from sub compJoinCDs ........ \n";

    print "\nThe size of the cd range array from sub compJoinCDs is = ", $length = @array_cdranges, "\n";

    $range_counter = 0;

    foreach $range (@array_cdranges ){
        $range_counter += 1;
#	my $cd_start;
#	my $cd_stop;

	my ($start, $stop) = split(/\.\./, $range);
	print "inside range loop in sub compJoinCDs, start position = ", $start, ", stop position = ", $stop, "\n";

	if($range_counter == 1){

	    $cd_start = $start;
	    $cd_stop = $stop;
	}
       
	my $start_at = $start -1; 
	
	print "in sub compJoinCDs: cds start - 1 due to index numbering = $start_at\n";
	
	print "in sub compJoinCDs, part of this cd fragment length is = ", $stop - $start_at, "\n";

	 $cds_sq .= substr ($cont_seq, $start_at, $stop - $start_at );

    }
    
    $cds_sq =~ s/''|NULL|\W//g;   
    $rev_cds_sq = reverse $cds_sq;    
    $rev_cds_sq =~ tr/ACGTacgt/TGCAtgca/;
	
#    print "rev comp seq in sub = ", $rev_cds_sq, "\n";
	
    print "After joing cds fragments from all ranges, the full cds in compJoinCds sub is:\n",$rev_cds_sq, "\n\n";
 
#    my $full_protein = translateProtein ($cds_sq);
#    print "\n\n****Printing from sub the full protein is:\n", $full_protein, "\n\n";       
    print "printing from compJoinCds sub  ended!!, \n";
 
 
    return $cd_start.":::".$rev_cds_sq.":::".$cd_stop;
    
 }

sub joinCds{

    my ($cd_ranges, $cont_seq) = @_;
    my $cds_sq;
    my $cd_start;
    my $cd_stop;
    
    ($cd_ranges_set  =  $cd_ranges) =~ s/\+|\*|join|forward/$1/g;
    
    my @array_cdranges = split (/\,/, $cd_ranges_set);

    print "printing from join CDs sub........ \n";

    print "\nThe size of the cd range array from sub joinCDs is = ", $length = @array_cdranges, "\n";

    $range_counter = 0;

    foreach $range (@array_cdranges ){
        $range_counter += 1;
#	my $cd_start;
#	my $cd_stop;

	my ($start, $stop) = split(/\.\./, $range);
	print "inside range loop in sub joinCDs, start position = ", $start, ", stop position = ", $stop, "\n";

	if($range_counter == 1){

	    $cd_start = $start;
	    $cd_stop = $stop;
	}
       
	my $start_at = $start -1; 
	
	print "in sub joinCDs: cds start - 1 due to index numbering = $start_at\n";
	
	print "in sub joinCDs, part of this cd fragment length is = ", $stop - $start_at, "\n";

	 $cds_sq .= substr ($cont_seq, $start_at, $stop - $start_at );

    }
    
    $cds_sq =~ s/''|NULL|\W//g;

    print "After joing cds fragments from all ranges, the full cds is:\n",$cds_sq, "\n\n";
 
#    my $full_protein = translateProtein ($cds_sq);
#    print "\n\n****Printing from sub the full protein is:\n", $full_protein, "\n\n";       
    print "printing from join CDs ended!!, \n";
 
   return $cd_start.":::".$cds_sq.":::".$cd_stop;
    
 }



sub cdsSeq{
    my ($cds, $cont_seq) = @_;
    my ($cds_sq);
    
    $cds =~ s/\s|\n|NULL//g;
    
#    print "cds range in sub for forward strand is = ", $cds, "\n";
    
    my ($cds_start, $cds_end) = split( /\.\./, $cds);
	
    $cds_start =~ s/\s|\n|NULL//g; 
    $cds_end =~ s/\s|\n|NULL//g;
    
    my  $cds_start_at = $cds_start -1;
    my  $cds_end_at = $cds_end -1;
    
#    $cds_sq = substr ($cont_seq, $cds_start_at, $cds_end_at - $cds_start_at );
    $cds_sq = substr ($cont_seq, $cds_start_at, $cds_end - $cds_start_at );

    print "printing from  sub cdsSeq for forward strand...... \n";
    print "start position = ", $cds_start, "\n";
    print "end position = ", $cds_end, "\n";
    
    print "cds seq length in cdsSeq sub is=", length $cds_sq, "\n";
    
    print "*** sequence of cds in cdsSeq sub is: ", $cds_sq, "\n";
    print"Finished printing on forward strand from sub cdsSeq....\n";

    return $cds_start.":::".$cds_sq.":::".$cds_end;

}

# ------------------------------------------------------------------------------------
# TABLES FOR METAGENOMICS DB.....................................................
# .....................................................................................

sub creatMetcontigTable{
    
    my($contig_name_ref, $species_ref, $taxonomy_ref, $form_ref, $contig_gc_ref, $contig_size_ref, $taxon_id_ref) = @_;
    $contig_name = $contig_name_ref;
    $species = $species_ref;
    $taxonomy = $taxonomy_ref;
    $form = $form_ref;
    $contig_gc = $contig_gc_ref;
    $contig_size = $contig_size_ref;
    $taxon_id = $taxon_id_ref;

#    $query_handle = $connect_to_host->prepare('INSERT into contigs (name, size, species, taxon_id, taxonomy, form, proteins, gc)VALUES (?,?,?,?,?,?,?,?)');
    $query_handle = $connect_to_host->prepare('INSERT into contigs (name, species, taxonomy, form, gc, size, taxon_id)VALUES (?,?,?,?,?,?,?)');

    $query_handle->execute($contig_name, $species, $taxonomy, $form, $contig_gc, $contig_size, $taxon_id) or die "could not execute statement: ".$query_handle;
    
    

}


sub creatMetProteinTable{#($pepseq, $contig_name);

    my ( $pepseq_ref, $contig_name_ref, $protein_name_ref) = @_;

    $protein_name = $protein_name_ref;
    $protein_sequence = $pepseq_ref;
    $contig_name = $contig_name_ref;


    $query_handle = $connect_to_host->prepare('SELECT id, name FROM contigs WHERE name = ?');
    $query_handle->execute($contig_name) or die "could not execute statement in sub createMetProteinTable: ".$query_handle;
    
    @data = $query_handle->fetchrow_array();
    
    print "The data extracted from contig's table in sub createProteinTable are:\n";
    print @data[0]."=".@data[1],"\n"; # @data[1]= name of contig, @data[0]= contig_id from contig table
    
    my $contig_id = @data[0];

    $query_handle = $connect_to_host->prepare('INSERT INTO proteins  (name, contig_id, seq)VALUES (?,?,? )');
    $query_handle->execute($protein_name,@data[0],$protein_sequence) or die "could not execute statement in sub createMetProteinTable: ".$query_handle;
    
    print " I have successfully inserted protein sequence info in the sub createMetProteinTable\n";
    
    
    $query_handle = $connect_to_host->prepare('SELECT id, name FROM proteins WHERE name = ?');
    $query_handle->execute($protein_name) or die "could not execute statement in sub createMetProteinTable: ".$query_handle;

    @data_protein = $query_handle->fetchrow_array();
    
    print "The protein data extracted from proteins table in sub createMetProteinTable are:\n";
    print "protein id = ", @data_protein[0]."="."protein name = ", @data_protein[1], "\n";
    
    return @data_protein[0];
}

sub creatMetCDSseqsTable {#($protein_id,$cdssequence);
    
    my ($protein_id_ref,$cds_seq_ref ) = @_;
    my $protein_id = $protein_id_ref;
    my $cds_seq = $cds_seq_ref;

    $query_handle = $connect_to_host->prepare('INSERT INTO  cdsseqs (protein_id, seq) VALUES (?,? )');
    $query_handle->execute($protein_id, $cds_seq ) or die "could not execute statementin sub createMetCDSseqTable : ".$query_handle;
    
    print "SUCESSFULLY ENTERED THE CEDS SEQ IN CDSSEQ TABLE\n";

}





