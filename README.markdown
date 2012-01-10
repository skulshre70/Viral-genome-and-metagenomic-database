>##README##
>
>A Software Program/Tool for Creating Viral and Metagenomic Database  
>=================================================================== 
>
>Author: Sanjeev Kulshreshtha (skulshre@yahoo.com)
>Date:   08/15/2011
>
>##WHY THIS PROGRAM?##
>
>There was a need to create a database at JCVI for phyologenetic analysis of  
>metagenomic sequences which was required to integrate with the umbrealla APIS  
>database at JCVI for comparative analysis.  
>  
>##WHAT THIS PROGRAM DOES?##
>
>This software program automatically creates a database comprising viral genomic  
>and microbial metagenomic expressed  sequence data derived from various public  
>biological databases.  The database can be viewed using SQL queries either using  
>a GUI interface or directly on the Unix, linux or PC kernels.  
>
>##HOW DOES THIS PROGRAM WORK?##
>
>1.  This program extracts, parses and generates final data that are loaded int  
>database.  The data are processed and validated with pre-defined criteria before  
>transporting to the database.  Only if the published peptide data and its own  
>translated peptide sequence data match, then only the data are loaded into  
>database.  
>
>2.  The database schema and format which is used by this program is MySQL.  
>This program uses Perl-DBI, Perl-SQL, Perl-MySQL, Perl and BioPerl modules.  
>The hall marks of the program are that-  
> >i) Validation of data enormously helps in exceptionally very high efficiency and  
>performance of the code in memory usage and speed.  
> >ii)  It avoids thread crowding and burdening database server.  
>
>3.  This program has a small main program, while majority of the work is modeled  
>using subroutines.  
>
>##WHY I AM PROUD OF IT?##
>
>The program works very fast and efficiently without slowing down grid network speed.  
