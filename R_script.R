#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$			REIDENTIFY			$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#A script for estimating the original identities between nucleic acid sequences (prior to introduction of sequencing errors)

##########################################################
#Specify file paths, user defined-variables, and packages#
##########################################################
#During normal use, change code under this heading only

#-------------------------------------------------
#Define paths for sequence (fasta or fastq) files)
#-------------------------------------------------
#All sequences must be trimmed to same region of interest (e.g., V4) for downstream alignment

#Set working directory
#Uncomment if not set by external script
#setwd("L:/REIDENTIFY/run_130401/data/")

#File path for sample sequences (those containing errors) 
#Uncomment if not set by external script
#sample.filepath=file.path=("Mock1_S1_L001_Merged_001_V34.fq")

#File path for reference sequences (those containing no errors)
#Uncomment if not set by external script
#reference.filepath=file.path=("reference_sequences_V34_no_primers.fasta") 

#File path for UPARSE-REF annotations (generated from USEARCH)
#Uncomment if not set by external script
#annotations.filepath=file.path=("Mock1_S1_L001_Merged_001_V34.up")

#------------------------------------------------------------------------
#Define maximum number of sequences to be analyzed and datapoints to plot
#------------------------------------------------------------------------
n_sample_max.scalar=1000 #Maximum number of sequences to be analyzed
n_datapoints_max.scalar=10000 #Maximum number of datapoints to plot
random_seed.binary=FALSE #For reproducibility, same sequences and datapoints are subsampled each time (set to TRUE if subsampling should be different)

#-------------------------
#Install and load packages
#-------------------------
#Install ShortRead package (uncomment below unless already installed)
#source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
#sasbiocLite("ShortRead")

#Install BioStrings package (uncomment below unless already installed)
#source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
#biocLite("Biostrings")

#Install msa package (uncomment below unless already installed)
#source("https://bioconductor.org/biocLite.R") #try http:// if https:// URLs are not supported
#biocLite("msa")

#Install svMisc package (uncomment below unless already installed)
#install.packages("svMisc")

#Install ggidentities_original_estimated_vs_original_actual.plot package (uncomment below unless already installed)
#install.packages("ggidentities_original_estimated_vs_original_actual.plot")

#Install gridExtra package (uncomment below unless already installed)
#install.packages("gridExtra")

#Install extrafont package and import fonts (uncomment below unless already installed)
#install.packages("extrafont")
#font_import()

#Load ShortRead package
library(ShortRead)

#Load Biostrings package
library(Biostrings)

#Load msa package
library(msa)

#Load svMisc package
library(svMisc)

#Load ggplot package
library(ggplot2)

#Load gridExtra package
library("gridExtra")

#Load extrafont package
library("extrafont")

#------------------------------------
#Load sequence (fasta or fastq) files
#------------------------------------

#Read in fastq file for sample sequences (not filtered for chimeras)
sample_unfiltered.shortreadq=readFastq(sample.filepath)

#Read in fasta file for reference sequences
reference.shortread=readFasta(reference.filepath)

#Read in fasta file for UPARSE-REF annotations
annotations.table=read.table(annotations.filepath, sep="\t", header=FALSE)

#-------------------------------------------------
#Set seed for subsampling sequences and datapoints
#-------------------------------------------------

if(random_seed.binary==TRUE) #Seed for random number generator itself is random
{
	rand=sample(0:10000,1) #Seed is between 0 to 10000
}
if(random_seed.binary==FALSE) #Seed for random number generator is same each time
{
	rand=1 #Seed is 1
}

###########################################################
#Remove sample sequences matching chimera sequences sample#
###########################################################

#Reorder UPARSE-REF annotations to match order of sample sequence IDs
annotations.dataframe=data.frame(annotations.table) #Place UPARSE-REF annotations in dataframe
target.vector=id(sample_unfiltered.shortreadq) #Place vsample sequence IDs in vector (target for match)
annotations_ordered.table=annotations.dataframe[match(target.vector, annotations.dataframe$V1),] #Match annotations in dataframe with IDs in target

#--------------------------------------------------------------------
#Use loop to determine which sample sequences match chimera sequences
#--------------------------------------------------------------------
#Define or initialize values of loop
n_sample_unfiltered.scalar=length(sample_unfiltered.shortreadq) #Number of sample sequences (not filtered for chimeras)
nochimera.vector=array(0,dim=c(0)) #Create empty array to hold index values of sequences that are not chimeras

#Perform loop
for(i in c(1:n_sample_unfiltered.scalar)) #Perform loop for i=1..n_sample_unfiltered.scalar
{
	if(annotations_ordered.table[i,2]=="perfect") #If uparse_ref annotation is that of a perfect non-chimeric sequence
	{
		#Update nochimera.vector
		nochimera.vector=cbind(nochimera.vector,i) #Sample sequence i added to nochimera.vector
	}

	if(annotations_ordered.table[i,2]=="good") #If uparse_ref annotation is that of a good non-chimeric sequence
	{
		#Update nochimera.vector
		nochimera.vector=cbind(nochimera.vector,i) #Sample sequence i added to nochimera.vector
	}

	if(annotations_ordered.table[i,2]=="noisy") #If uparse_ref annotation is that of a noisy non-chimeric sequence
	{
		#Update nochimera.vector
		nochimera.vector=cbind(nochimera.vector,i) #Sample sequence i added to nochimera.vector
	}

	if(annotations_ordered.table[i,2]=="other") #If uparse_ref annotation is that of a other non-chimeric sequence
	{
		#Update nochimera.vector
		nochimera.vector=cbind(nochimera.vector,i) #Sample sequence i added to nochimera.vector
	}

	if(annotations_ordered.table[i,2]=="perfect_chimera")	#If uparse_ref annotation is that of a perfect chimera sequence
	{
		#Update nochimera.vector
		nochimera.vector=nochimera.vector #Sample sequence not added to nochimera.vector (nochimera.vector not updated)
	}
	
	if(annotations_ordered.table[i,2]=="noisy_chimera") #If uparse_ref annotation is that of a noisy chimera sequence
	{
		#Update nochimera.vector
		nochimera.vector=nochimera.vector #Sample sequence not added to nochimera.vector (nochimera.vector not updated)
	}

#Show progress of loop
progress(value=(i-1),max.value=n_sample_unfiltered.scalar) 
}

#-----------------------------------------------------------------
#Retain only those sample sequences not matching chimera sequences
#-----------------------------------------------------------------
sample_filtered.shortreadq=sample_unfiltered.shortreadq[as.vector(nochimera.vector)]

#################################################################
#Remove sample sequences if number exceeds user-defined maximum #
#################################################################
#Takes a random subsample sample sequences before removal

n_sample_filtered.scalar=length(sample_filtered.shortreadq) #Determine number of sequences after filtering for chimeras

if(n_sample_filtered.scalar>n_sample_max.scalar) #If number of sequences is more than the maximum
{
	set.seed(rand);random_number_sequence.vector=sample(1:n_sample_filtered.scalar, size=n_sample_max.scalar) #Create a vector of n_sample_max.scalar size with elements being randomly sampled from n_sample_filtered.scalar
	
	sample.shortreadq=sample_filtered.shortreadq[random_number_sequence.vector] #Include only first n_sample_max.scalar sequences 
}
if(n_sample_filtered.scalar<=n_sample_max.scalar) #If number of sequences is less or equal to the maximum
{
	sample.shortreadq=sample_filtered.shortreadq #Retain all sequences
}

#############################################################################
#Determine which reference sequence has highest identities to sample sequences#
#############################################################################

#--------------------------------------------------------
#Align sample and reference sequences (with ClustalOmega)
#--------------------------------------------------------
#Define number of sequences
n_sample.scalar=length(sample.shortreadq) #Number of sample sequences
n_reference.scalar=length(reference.shortread) #Number of reference sequences

#Combine sample and reference sequences into one set
combined.dnastring=xscat(c(sread(sample.shortreadq),sread(reference.shortread)))

#Perform alignment
alignment.msa=msa(combined.dnastring, method="ClustalOmega", type="dna", order="input")

#Output alignment as matrix
alignment.matrix=as.matrix(alignment.msa)

#Split matrix into one containing sample sequences and another containing reference sequences 
letters_sample_aligned.matrix=alignment.matrix[0:n_sample.scalar,]
letters_reference_aligned.matrix=alignment.matrix[(n_sample.scalar+1):(n_sample.scalar+n_reference.scalar),]

#-----------------------------------------------------------------------------------------------------------------
#Use loop to determine reference sequence with highest identities to sample sequence (and value of highest identities)
#-----------------------------------------------------------------------------------------------------------------
#Define or initialize values of loop
n_align.scalar=ncol(alignment.matrix) #The width of the alignment (number of nucleotide positions in the aligned sequences)
identities_referencematch.matrix=array(0,dim=c(n_sample.scalar))#Create empty array to hold identities values, with n_sample.scalar dimensions 
reference_letters.vector=array(0,dim=c(n_sample.scalar,n_align.scalar))#Create empty array to hold letters of reference sequence that align to the sample sequence, with n_sample.scalar x n_align.scalar dimensions

for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
{

	for(j in c(1:n_reference.scalar)) #Loop for j=1..n_reference.scalar
	{
		#*****************************************************************
		#Perform subloop to determine identities at each nucleotide position
		#*****************************************************************
		#Initialize values for subloop
		identities_referencematch.vector=array(0,dim=n_align.scalar) #Creates blank array to hold observed identities values at each nucleotide position, with n_align.scalar dimensions
		n_indel.scalar=0 #Number of indels shared between sample and reference sequences at a given k

		for(k in c(1:n_align.scalar)) #Loop for k=1..n_align.scalar
		{
			if((letters_sample_aligned.matrix[i,k]!="-")&(letters_reference_aligned.matrix[j,k]!="-")) #Calculate identities if no indels for letters_sample_aligned.matrix[i,k] and letters_reference_aligned.matrix[i,k]
			{
				#Calculate identities
				identities_referencematch.vector[k]=letters_sample_aligned.matrix[i,k]==letters_reference_aligned.matrix[j,k] 
			}
			else if((letters_sample_aligned.matrix[i,k]=="-")&(letters_reference_aligned.matrix[j,k]=="-")) #Set identities to 0 if indels for both letters_sample_aligned.matrix[i,k] and letters_reference_aligned.matrix[i,k]
			{
				#Calculate identities
				identities_referencematch.vector[k]=0
				
				#Count shared indels
				n_indel.scalar=n_indel.scalar+1 
			}
			else if(letters_sample_aligned.matrix[i,k]=="-")  #Set identities to 0 if indels at letters_sample_aligned.matrix[i,k] 		
			{
				#Calculate identities
				identities_referencematch.vector[k]=0 
			}
			else if(letters_reference_aligned.matrix[j,k]=="-")  #Set identities to 0 if indels at letters_reference_aligned.matrix[j,k] alone		
			{
				#Calculate identities
				identities_referencematch.vector[k]=0 
			}
		}

		identities_referencematch_max.scalar=sum(identities_referencematch.vector)/(n_align.scalar-n_indel.scalar) #Calculate observed identities for i and j over k
     
		#******************************************
		#Select for alignment with highest identities
		#******************************************
		if(identities_referencematch_max.scalar>identities_referencematch.matrix[i]) 
		{
			identities_referencematch.matrix[i]=identities_referencematch_max.scalar #Update value of highest identities 
			reference_letters.vector[i,]=letters_reference_aligned.matrix[j,] #Output letters of reference sequence with highest identities to sample sequence
		}
	
		#Show progress of loop
		progress(value=j+(i-1)*n_reference.scalar,max.value=n_sample.scalar*n_reference.scalar) 		

		#Advance indices
		if(identities_referencematch.matrix[i]==1) #Break loop if reference sequence with perfect identities is found
		{
			break
		}
	}
}

###################################################
#Calculate error probabilities from quality scores# 
###################################################

#-------------------------------------------------------------------------------------------
#Use loop to calculate error probabilities and bring into same alignment as sequence letters
#-------------------------------------------------------------------------------------------
#Defined or initialize values for loop
n_align.scalar=ncol(alignment.matrix) #The width of the alignment(number of nucleotide positions in the aligned sequences)
errorrates_sample_aligned.matrix=array(0,dim=c(n_sample.scalar, n_align.scalar)) #Create empty array to hold aligned error probabilities, with n_sample.scalar x n_align.scalar dimensions
i=1 #Index for error rates of sequence i

#Perform loop
for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
{
	#Calculate error rates
	errorrates_sample.matrix=10^(-(as(quality(sample.shortreadq)[i],"matrix"))/10) 

	#***************************************************************************
	#Perform subloop to determine error probabilites at each nucleotide position
	#***************************************************************************
	#Initialize values for subloop
	k=1 #Index for nucleotide position in aligned error probabilities
	l=1 #Index for nucleotide position in unaligned error probabilities

	for(k in c(1:n_align.scalar)) #Perform loop for k=1..n_align.scalar 
	{
		if(letters_sample_aligned.matrix[i,k]!="-") #Set error probability to errorrates_sample.matrix[l] if no indels
			{
				errorrates_sample_aligned.matrix[i,k]=errorrates_sample.matrix[l] 
			
				#Advance indices
				l=l+1
			}

		else if(letters_sample_aligned.matrix[i,k]=="-") #Set error probability to 0 if indels for letters_sample_aligned.matrix[i,k]
			{
				errorrates_sample_aligned.matrix[i,k]=0 
			
				#Advance indices
				l=l #Do not advance l because of indel
			}
	}
}

#############################################################################################
#Calculate observed and orignal identities (estimated) between all pairs of sample sequences#
#############################################################################################

#----------------------------------
#Create matrices to hold identities 
#----------------------------------
#Create matrices with n_sample.scalar diagonal elements set to 1
identities_observed.matrix=diag(n_sample.scalar) #Matrix for observed identities values
identities_original.matrix=diag(n_sample.scalar) #Create matrix for original identities values

#------------------------------------------------------
#Use loop to determine identities above matrix diagonal
#------------------------------------------------------
for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
{
	if(i<n_sample.scalar)
	{
		for(j in c((i+1):n_sample.scalar)) #Perform loop for j=1..(i-1)
		{
			#*****************************************************************
			#Perform subloop to determine identities at each nucleotide position
			#*****************************************************************
			#Initialize values for subloop
			identities_observed.vector=array(0,dim=n_align.scalar) #Creates blank array to hold observed identities values at each nucleotide position, with n_align.scalar dimensions
			identities_original_estimated.vector=array(0,dim=n_align.scalar) #Creates blank array to hold original identities values at each nucleotide position, with n_align.scalar dimensions
			n_indel.scalar=0 #Number of indels shared between sample and reference sequences at a given k

			for(k in c(1:n_align.scalar)) #Perform subloop for k=1..n_align.scalar
			{
				if((letters_sample_aligned.matrix[i,k]!="-")&(letters_sample_aligned.matrix[j,k]!="-")) #Calculate identities if no indels for letters_sample_aligned.matrix[i,k] and letters_sample_aligned.matrix[j,k]
				{
					#Calculate identities
					identities_observed.vector[k]=letters_sample_aligned.matrix[i,k]==letters_sample_aligned.matrix[j,k] #Observed identities
					identities_original_estimated.vector[k]=(9*identities_observed.vector[k]-3*errorrates_sample_aligned.matrix[i,k]-3*errorrates_sample_aligned.matrix[j,k]+4*(errorrates_sample_aligned.matrix[i,k]*errorrates_sample_aligned.matrix[j,k]))/(9-12*errorrates_sample_aligned.matrix[i,k]-12*errorrates_sample_aligned.matrix[j,k]+16*(errorrates_sample_aligned.matrix[i,k]*errorrates_sample_aligned.matrix[j,k])) #Original identities (see equation in manuscript)
				}
				else if((letters_sample_aligned.matrix[i,k]=="-")&(letters_sample_aligned.matrix[j,k]=="-")) #Set identities to 0 if indels for both letters_sample_aligned.matrix[i,k] and letters_sample_aligned.matrix[j,k]
				{
					#Calculate identities
					identities_observed.vector[k]=0 #Observed identities
					identities_original_estimated.vector[k]=0 #Original identities (see equation in manuscript)
		
					#Count shared indels
					n_indel.scalar=n_indel.scalar+1 
				}
				else if(letters_sample_aligned.matrix[i,k]=="-")  #Set identities to 0 if indels at letters_sample_aligned.matrix[i,k] alone, but do not increase count of shared indels		
				{
					#Calculate identities
					identities_observed.vector[k]=0 #Observed identities
					identities_original_estimated.vector[k]=0 #Original identities

					#Count shared indels
					n_indel.scalar=n_indel.scalar				
				}
				else if(letters_sample_aligned.matrix[j,k]=="-")  #Set identities to 0 if indels at letters_sample_aligned.matrix[j,k] alone, but do not increase count of shared indels		
				{
					#Calculate identities
					identities_observed.vector[k]=0 #Observed identities
					identities_original_estimated.vector[k]=0 #Original identities

					#Count shared indels
					n_indel.scalar=n_indel.scalar				
				}
			}
			#Show progress of loop
			progress(value=j+(i-1)*n_sample.scalar,max.value=n_sample.scalar*(n_sample.scalar-1)) 

			#Calculates identities across nucleotide positions for sequences i and j
			identities_observed.matrix[i,j]=sum(identities_observed.vector)/(n_align.scalar-n_indel.scalar) #Observed identities
			identities_original.matrix[i,j]=sum(identities_original_estimated.vector)/(n_align.scalar-n_indel.scalar) #Original identities
		}
	}
}

#----------------------------------------
#Fill in identities below matrix diagonal
#----------------------------------------
identities_observed.matrix[lower.tri(identities_observed.matrix)] = t(identities_observed.matrix)[lower.tri(identities_observed.matrix)] #Observed identities
identities_original.matrix[lower.tri(identities_original.matrix)] = t(identities_original.matrix)[lower.tri(identities_original.matrix)] #Original identities

#####################################################################################################
#Determine the original identities (actual) between references sequences that match sample sequences#
#####################################################################################################

#----------------------------------
#Create matrices to hold identities 
#----------------------------------
#Create matrices with n_sample.scalar diagonal elements set to 1
identities_original_actual.matrix=diag(n_sample.scalar) #Matrix for identities values

#------------------------------------------------------
#Use loop to determine identities above matrix diagonal
#------------------------------------------------------
for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
{
	if(i<n_sample.scalar)
	{
		for(j in c((i+1):n_sample.scalar)) #Perform loop for j=1..(i-1)
		{
			#*****************************************************************
			#Perform subloop to determine identities at each nucleotide position
			#*****************************************************************
			#Initialize values for subloop
			identities_original_actual.vector=array(0,dim=n_align.scalar) #Creates blank array to hold observed identities values at each nucleotide position, with n_align.scalar dimensions
			n_indel.scalar=0 #Number of indels shared between sample and reference sequences at a given k

			#Perform subloop
			for(k in c(1:n_align.scalar)) #Perform subloop for k=1..n_align.scalar
			{
				if((reference_letters.vector[i,k]!="-")&(reference_letters.vector[j,k]!="-")) #Calculate identities if no indels for reference_letters.vector[i,k] and reference_letters.vector[j,k]
				{
					#Calculate identities
					identities_original_actual.vector[k]=reference_letters.vector[i,k]==reference_letters.vector[j,k] #Observed identities
				}
				else if((reference_letters.vector[i,k]=="-")&(reference_letters.vector[j,k]=="-")) #Set identities to 0 if indels for both reference_letters.vector[i,k] and reference_letters.vector[j,k]
				{
					#Calculate identities
					identities_original_actual.vector[k]=0 #Observed identities
		
					#Count shared indels
					n_indel.scalar=n_indel.scalar+1 
				}
				else if(reference_letters.vector[i,k]=="-")  #Set identities to 0 if indels at reference_letters.vector[i,k] alone, but do not increase count of shared indels		
				{
					#Calculate identities
					identities_original_actual.vector[k]=0 #Observed identities

					#Count shared indels
					n_indel.scalar=n_indel.scalar				
				}
				else if(reference_letters.vector[j,k]=="-")  #Set identities to 0 if indels at reference_letters.vector[j,k] alone, but do not increase count of shared indels		
				{
					#Calculate identities
					identities_original_actual.vector[k]=0 #Observed identities

					#Count shared indels
					n_indel.scalar=n_indel.scalar					
				}
			}
			#Show progress of loop
			progress(value=j+(i-1)*n_sample.scalar,max.value=n_sample.scalar*n_sample.scalar) 

			#Calculates identities across nucleotide positions for sequences i and j
			identities_original_actual.matrix[i,j]=sum(identities_original_actual.vector)/(n_align.scalar-n_indel.scalar) 
		}
	}
}
#----------------------------------------
#Fill in identities below matrix diagonal
#----------------------------------------
identities_original_actual.matrix[lower.tri(identities_original_actual.matrix)] = t(identities_original_actual.matrix)[lower.tri(identities_original_actual.matrix)] #Observed identities

#######################################################################################
#Calculate actual sequencing error rate and compare to instrument-reported error rates#
#######################################################################################

#-------------------------------------------------------------------------
#Calculate identities between sample sequence and reference sequence match
#-------------------------------------------------------------------------
#Define or initialize values for loop
errorrates_actual.vector=array(0,dim=c(n_sample.scalar,n_align.scalar)) #Creates empty array to hold error rates at each nucleotide position, with n_sample.scalar x n_align.scalar dimensions

for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
{
		#*****************************************************************
		#Perform subloop to determine errors at each nucleotide position
		#*****************************************************************
		#Initialize values for subloop

		#Perform subloop
		for(k in c(1:n_align.scalar)) #Perform subloop for k=1..n_align.scalar
		{
			if((letters_sample_aligned.matrix[i,k]!="-")&(reference_letters.vector[i,k]!="-")) #Calculate error rate if no indels for letters_sample_aligned.matrix[i,k] and reference_letters.vector[i,k]
			{
				#Calculate identities
				errorrates_actual.vector[i,k]=1-(letters_sample_aligned.matrix[i,k]==reference_letters.vector[i,k]) #Set rate rate to 0 if no match and 1 if match
			}
			else if((letters_sample_aligned.matrix[i,k]=="-")&(reference_letters.vector[i,k]=="-")) #Set error rate to NA if indels for both letters_sample_aligned.matrix[i,k] and reference_letters.vector[i,k]
			{
				#Calculate identities
				errorrates_actual.vector[i,k]=NA 
			}
			else if(letters_sample_aligned.matrix[i,k]=="-")  #Set error rate to 1 if indels at letters_sample_aligned.matrix[i,k] alone	
			{
				#Calculate identities
				errorrates_actual.vector[i,k]=1 
			}
			else if(reference_letters.vector[i,k]=="-")  #Set error rate to 1 if indels at reference_letters.vector[i,k] alone	
			{
				#Calculate identities
				errorrates_actual.vector[i,k]=1 			
			}
		}
		#Show progress of loop
		progress(value=(i-1)*n_sample.scalar,max.value=n_sample.scalar*n_sample.scalar) 
}

#------------------------------------------------------
#Output actual and predicted error rates into dataframe
#------------------------------------------------------

#Output errorrates_actual.vector[i,k] and errorrates_sample_aligned.matrix[i,k] into columns of 2-dimensional array
errorrates.dataframe=array(NA,dim=c(n_sample.scalar*n_align.scalar,2))

for(i in c(1:n_sample.scalar)) #Perform loop for i=1..n_sample.scalar
{
	for(k in c(1:n_align.scalar)) #Perform subloop for k=1..n_align.scalar
	{
		errorrates.dataframe[(k+(i-1)*n_align.scalar),1]=errorrates_actual.vector[i,k] #
		errorrates.dataframe[(k+(i-1)*n_align.scalar),2]=errorrates_sample_aligned.matrix[i,k] #
	}
}

#Remove rows containing "NA"
errorrates.matrix=as.matrix(na.omit(errorrates.dataframe))

###############################################################################
#Calculate RMSPE, bias, actual error rates, and instrument-reported error rates
###############################################################################
mspe_identities_observed.scalar=(mean((identities_original_actual.matrix-identities_observed.matrix)^2))^0.5
bias_identities_observed.scalar=mean(identities_original_actual.matrix-identities_observed.matrix)

mspe_identities_original.scalar=(mean((identities_original_actual.matrix-identities_original.matrix)^2))^0.5
bias_identities_original.scalar=mean(identities_original_actual.matrix-identities_original.matrix)

errorrates_actual.scalar=mean(errorrates.matrix[,1])
errorrates_sample.scalar=mean(errorrates.matrix[,2])

##############
#Create plots#
##############

#Load fonts for plots
loadfonts(device = "win")

#--------------------------------------------
#Plot observed vs. actual observed identities
#--------------------------------------------

#Output identities_observed.matrix and identities_original.matrix into rows of a dataframe  
identities_observed.vector=as.vector(identities_observed.matrix) #Output identities_observed.matrix into a vector
identities_observed_nodiagonals.vector=identities_observed.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
identities_original_actual.vector=as.vector(identities_original_actual.matrix) #Output identities_original_actual.matrix into a vector
identities_original_actual_nodiagonals.vector=identities_original_actual.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
identities_observed_vs_original_actual.matrix=cbind(identities_observed_nodiagonals.vector,identities_original_actual_nodiagonals.vector)
identities_observed_vs_original_actual.dataframe=data.frame(identities_observed_vs_original_actual.matrix)

#**********************************************************************
#Remove datapoints from being plotted if exceeding user-defined maximum
#**********************************************************************

n_datapoints.scalar=length(sample.shortreadq)*(length(sample.shortreadq)-1) #Number of datapoints

if(n_datapoints.scalar>n_datapoints_max.scalar) #If number of datapoints is more than the maximum
{
	set.seed(rand);random_number_datapoint.vector=sample((1:n_datapoints.scalar), size=n_datapoints_max.scalar) #Create a vector of n_datapoints_max.scalar size with elements being randomly sampled from n_datapoints.scalar

	identities_observed_vs_original_actual.dataframe=identities_observed_vs_original_actual.dataframe[random_number_datapoint.vector,] #Include only first n_datapoints.scalar datapoints 	
}
if(n_datapoints.scalar<=n_datapoints_max.scalar) #If number of datapoints is less or equal to the maximum
{
	identities_observed_vs_original_actual.dataframe=identities_observed_vs_original_actual.dataframe #Retain all datapoints
}

#********************************************
#Make dataframe for plotting 1:1 (unity) line
#********************************************
unity_line.dataframe=data.frame(x=c(60,105),y=c(60,105))

#*********
#Make plot
#*********

identities_observed_vs_original_actual.plot=ggplot(data=identities_observed_vs_original_actual.dataframe,aes(x=identities_original_actual_nodiagonals.vector*100,y=identities_observed_nodiagonals.vector*100)) #Specify dataframe and x and y variables
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+geom_line(data=unity_line.dataframe,aes(x=x,y=y,color="1:1 line",size="1:1 line",linetype="1:1 line")) #Plot 1:1 line
marker_border.string=rgb(red=0, green=112, blue=192, maxColorValue = 255)#Specify marker border color
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+geom_point(col=marker_border.string, shape=1, alpha=0.1, size=0.25) #Set marker borders
marker_fill.string=rgb(red=0, green=112, blue=192, maxColorValue = 255)#Specify marker fill color
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+geom_point(col=marker_fill.string, shape=1, alpha=0.05, size=0.25) #Plot markers
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+geom_smooth(aes(color="Best fit line",size="Best fit line",linetype="Best fit line"),method="lm",se=FALSE) #Plot best-fit line (linear)
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+scale_color_manual(name="", values=c("grey80", "black")) #Specify line color
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+scale_size_manual(name="",values=c(0.5,0.25)) #Specify line size
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+scale_linetype_manual(name="",values=c("solid","88")) #Specify line type
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+labs(x="Original sequence identity (actual) (%)", y="Observed sequence identity (%)") #Define axis labels
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+xlim(60,105)+ylim(60,105) #Set x and y scale limits
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(fill = "transparent",colour = NA), legend.key= element_rect(fill = "transparent",colour = NA)) #Set panel, plot, legend, and legend key background to transparent
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+theme(axis.title=element_text(family="Arial", face="bold", colour="black", size="10"),axis.text=element_text(family="Arial",face="bold", color="black", size="10"),legend.text=element_text(family="Arial",face="bold", color="black", size="10")) #Define text settings
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+theme(axis.ticks=element_line(colour="gray80",size=0.5), axis.ticks.length=unit(-0.25, "cm"), axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) ) #Define axis tick mark settings
identities_observed_vs_original_actual.plot=identities_observed_vs_original_actual.plot+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray80", size=0.5)) #Remove gridlines and border; define axis line settings

#---------------------------------------------
#Plot estimated vs. actual observed identities
#---------------------------------------------

#Output identities_observed.matrix and identities_original.matrix into rows of a dataframe  

identities_original_estimated.vector=as.vector(identities_original.matrix) #Output identities_original.matrix into a vector
identities_original_estimated_nodiagonals.vector=identities_original_estimated.vector[-seq(1,n_sample.scalar^2,n_sample.scalar+1)] #Remove diagonal elements
identities_original_estimated_vs_original_actual.matrix=cbind(identities_original_estimated_nodiagonals.vector,identities_original_actual_nodiagonals.vector)
identities_original_estimated_vs_original_actual.dataframe=data.frame(identities_original_estimated_vs_original_actual.matrix)

#**********************************************************************
#Remove datapoints from being plotted if exceeding user-defined maximum
#**********************************************************************

n_datapoints.scalar=length(sample.shortreadq)*(length(sample.shortreadq)-1) #Number of datapoints

if(n_datapoints.scalar>n_datapoints_max.scalar) #If number of datapoints is more than the maximum
{
	random_number_datapoint.vector=sample((1:n_datapoints.scalar), size=n_datapoints_max.scalar) #Create a vector of n_datapoints_max.scalar size with elements being randomly sampled from n_datapoints.scalar

	identities_original_estimated_vs_original_actual.dataframe=identities_original_estimated_vs_original_actual.dataframe[random_number_datapoint.vector,] #Include only first n_datapoints.scalar datapoints 	
}
if(n_datapoints.scalar<=n_datapoints_max.scalar) #If number of datapoints is less or equal to the maximum
{
	identities_original_estimated_vs_original_actual.dataframe=identities_original_estimated_vs_original_actual.dataframe #Retain all datapoints
}

#*********
#Make plot
#*********

identities_original_estimated_vs_original_actual.plot=ggplot(data=identities_original_estimated_vs_original_actual.dataframe,aes(x=identities_original_actual_nodiagonals.vector*100, y=identities_original_estimated_nodiagonals.vector*100)) #Specify dataframe and x and y variables
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+geom_line(data=unity_line.dataframe,aes(x=x,y=y,color="1:1 line",size="1:1 line",linetype="1:1 line")) #Plot 1:1 line
marker_border.string=rgb(red=0, green=112, blue=192, maxColorValue = 255)#Specify marker border color
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+geom_point(col=marker_border.string, shape=1, alpha=0.1, size=0.25) #Set marker borders
marker_fill.string=rgb(red=0, green=112, blue=192, maxColorValue = 255)#Specify marker fill color
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+geom_point(col=marker_fill.string, shape=1, alpha=0.05, size=0.25) #Plot markers
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+geom_smooth(aes(color="Best fit line",size="Best fit line",linetype="Best fit line"),method="lm",se=FALSE) #Plot best-fit line (linear)
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+scale_color_manual(name="", values=c("grey80", "black")) #Specify line color
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+scale_size_manual(name="",values=c(0.5,0.25)) #Specify line size
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+scale_linetype_manual(name="",values=c("solid","88")) #Specify line type
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+labs(x="Original sequence identity (actual) (%)", y="Original sequence identity (estimated) (%)") #Define axis labels
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+xlim(60,105)+ylim(60,105) #Set x and y scale limits
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(fill = "transparent",colour = NA), legend.key= element_rect(fill = "transparent",colour = NA)) #Set panel, plot, legend, and legend key background to transparent
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+theme(axis.title=element_text(family="Arial", face="bold", colour="black", size="10"),axis.text=element_text(family="Arial",face="bold", color="black", size="10"),legend.text=element_text(family="Arial",face="bold", color="black", size="10")) #Define text settings
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+theme(axis.title=element_text(margin=unit(c(1,1,1,1), "cm"))) #Define margin around axis titles
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+theme(axis.ticks=element_line(colour="gray80",size=0.5), axis.ticks.length=unit(-0.25, "cm"), axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) ) #Define axis tick mark settings
identities_original_estimated_vs_original_actual.plot=identities_original_estimated_vs_original_actual.plot+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray80", size=0.5)) #Remove gridlines and border; define axis line settings

grid.arrange(identities_observed_vs_original_actual.plot, identities_original_estimated_vs_original_actual.plot, ncol=2)

#######################
#Output plots to *.pdf#
#######################
#Uncomment if plots should be outputted
#ggsave("identities_observed_vs_original_actual.pdf", device="pdf", plot=identities_observed_vs_original_actual.plot, width=5, height=3.5, units="in", bg="transparent") #Save file as pdf with dimensions and units as shown, as well as transparent background
#ggsave("identities_original_estimated_vs_original_actual.pdf", device="pdf", plot=identities_original_estimated_vs_original_actual.plot, width=5, height=3.5, units="in", bg="transparent") #Save file as pdf with dimensions and units as shown, as well as transparent background

#########################
#Output results to *.csv#
#########################
#Uncomment if results should be outputted

#Reference sequence with highest identities to sample sequence and its identities
#write.csv(reference_letters.vector, file = "reference_with_highest_identities_to_sample.csv") #Reference sequence
#write.csv(identities_referencematch.matrix, file = "identities_between_sample_and_reference.csv") #Identity

#Observed and original identities between pairs of sample sequences
#write.csv(identities_observed.matrix, file = "observed_identities.csv") #Observed identities
#write.csv(identities_original.matrix, file = "original_identities.csv") #Original identities

#Original identities between pairs of reference sequences matches
#write.csv(identities_original_actual.matrix, file = "observed_identities_reference_sequence_matches.csv")

#Observed, estimated original, and actual original identities (one per column)
#write.csv(cbind(identities_observed_nodiagonals.vector,identities_original_estimated_nodiagonals.vector,identities_original_actual_nodiagonals.vector), file = "identities_comparison.csv")
