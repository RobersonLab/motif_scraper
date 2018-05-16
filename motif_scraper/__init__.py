#!/usr/bin/env python

"""
motif_scraper: a tool for searching for nucleic acid motifs using degenerate nucleotide queries
"""

##########
# import #
##########
import multiprocessing as mp
import pyfaidx
import regex
import six
import sys
import time

####################
# Version and name #
####################
__script_path__ = sys.argv[0]
__script_name__ = __script_path__.split('/')[-1].split('\\')[-1]
__version__ = '1.0.1'

########
# fxns #
########	
def rev_comp( seq, molecule='dna' ):
	""" DNA|RNA seq -> reverse complement
	"""
	
	if molecule == 'dna':
		nuc_dict = { "A":"T", "B":"V", "C":"G", "D":"H", "G":"C", "H":"D", "K":"M", "M":"K", "N":"N", "R":"Y", "S":"S", "T":"A", "V":"B", "W":"W", "Y":"R" }
	elif molecule == 'rna':
		nuc_dict = { "A":"U", "B":"V", "C":"G", "D":"H", "G":"C", "H":"D", "K":"M", "M":"K", "N":"N", "R":"Y", "S":"S", "U":"A", "V":"B", "W":"W", "Y":"R" }
	else:
		raise ValueError( "rev_comp requires molecule to be dna or rna" )
	
	if not isinstance( seq, six.string_types ):
		raise TypeError( "seq must be a string!" )
	
	return ''.join( [ nuc_dict[c] for c in seq.upper()[::-1] ] )

def make_degenerate_regex( motif_seq, molecule='dna' ):
	""" Degenerate sequence -> regex
	Example: NNYCGAARN -> [ACGT]{2}[CT]CGA{2}[AG][ACGT]
	"""
	
	if not isinstance( motif_seq, six.string_types ):
		raise TypeError( "motif_seq must be a string!" )
		
	if molecule == 'dna':
		degenerate_code = { "A":"A", "B":"[CGT]", "C":"C", "D":"[AGT]", "G":"G", "H":"[ACT]", "K":"[GT]", "M":"[AC]", "N":"[ACGT]", "R":"[AG]", "S":"[GC]", "T":"T", "V":"[ACG]", "W":"[AT]", "Y":"[CT]" }
	elif molecule == 'rna':
		degenerate_code = { "A":"A", "B":"[CGU]", "C":"C", "D":"[AGU]", "G":"G", "H":"[ACU]", "K":"[GU]", "M":"[AC]", "N":"[ACGU]", "R":"[AG]", "S":"[GC]", "U":"U", "V":"[ACG]", "W":"[AU]", "Y":"[CU]" }
	else:
		raise ValueError( "make_degenerate_regex requires molecule to be dna or rna" )
	
	regex_string = ''
	
	idx = 0
	
	while idx < len( motif_seq ):
		curr = motif_seq[idx]
		
		count = 1
		
		for next_idx in range( idx+1, len( motif_seq ) ):
			next = motif_seq[next_idx]
			
			if next == curr:
				count += 1
			else:
				break
		
		regex_string += degenerate_code[curr]
		
		if count > 1:
			idx = idx + count - 1
			regex_string += "{%s}" % ( count )
		
		idx += 1
	
	return regex_string
	
###########
# classes #
###########
class SequenceMotif:
	"""
	motif = the search motif
	sequence = the matching nucleic acid sequence
	contig = name of the FASTA contig where match was found
	positionStart = first base of the match*
	strand = +|- strand of match relative to FASTA sequence
	regexMatch = the regex match
	
	*positionStart is expected to be 0-base. Very important!
	"""
	
	def __init__( self, motif, contig, positionStart, strand, regexMatch, molecule='dna', regex_match_start=None, regex_end_val=None, regex_group_sequence=None ):
		if regexMatch is not None:
			if strand == '+':
				self.seq = regexMatch.group()
				self.start = positionStart + regexMatch.start() + 1
				self.end = self.start + len( self.seq ) - 1
			elif strand == '-':
				self.seq = rev_comp( regexMatch.group(), molecule )
				self.end = positionStart + regexMatch.start() + 1
				self.start = self.end + len( self.seq ) - 1
			else:
				raise ValueError( "SequenceMotif classes require strand to be either '+' or '-'" )
		else:
			if strand == '+':
				self.seq = regex_group_sequence
				self.start = positionStart + regex_match_start
				self.end = self.start + len( self.seq ) - 1
			elif strand == '-':
				self.seq = rev_comp( regex_group_sequence, molecule )
				self.end = positionStart + regex_match_start
				self.start = self.end + len( self.seq ) - 1
			else:
				raise ValueError( "SequenceMotif classes require strand to be either '+' or '-'" )
		
		self.motif = motif
		self.contig = contig
		self.strand = strand
	
	def __str__( self ):
		""" string interpolation for file writing as a CSV
		"""
		
		return "%s,%s,%s,%s,%s,%s" % ( self.contig, self.start, self.end, self.strand, self.seq, self.motif )
		   
	def __repr__( self ):
		return str( self )
	
	def __eq__( self, other ):
		if self.motif == other.motif and self.contig == other.contig and self.strand == other.strand and self.start == other.start and self.end == other.end and self.seq == other.seq:
			return True
		else:
			return False

######################
# fxns using classes #
######################
def fasta_motif_scan( fasta_fname, input_tuples, regex_ready=False, allow_overlaps=True, molecule='dna'  ):
	"""
	fasta_fname = string path to FASTA file
	input_tuples = tuple containing (1) motif sequence, (2) contig name, (3) start position*, (4) end position, (5) strand to search
	
	*start is expected to be 0-base, end is expected to be 1-base
	"""
	
	###################
	# validity checks #
	###################
	if not isinstance( fasta_fname, six.string_types ):
		raise TypeError( "In fasta_motif_scan, fasta_fname must be a string!" )
	elif not isinstance( molecule, six.string_types ):
		raise TypeError( "In fasta_motif_scan, molecule must be a string!" )
	elif isinstance( input_tuples, six.string_types ) or not isinstance( input_tuples, tuple ):
		raise( TypeError( "In fasta_motif_scan, input_tuples should be a tuple!" ) )
	elif not type( regex_ready ) is bool:
		raise( TypeError( "In fasta_motif_scan, regex_ready should be a bool!" ) )
	
	#########################
	# setup some local vars #
	#########################
	motif_seq, contig, start, end, strand = input_tuples
	
	if regex_ready == False:
		if strand == '+':
			regex_compiled = regex.compile( make_degenerate_regex( motif_seq ) )
		elif strand == '-':
			regex_compiled = regex.compile( make_degenerate_regex( rev_comp( motif_seq, molecule ) ) )
	else:
		if strand == '+':
			regex_compiled = regex.compile( motif_seq )
		elif strand == '-':
			regex_compiled = regex.compile( rev_comp( motif_seq, molecule ) )
	
	site_list = []
	
	with pyfaidx.Fasta( fasta_fname, as_raw=True ) as FAIDX:
		sequence = str( FAIDX[contig][start:end] ).upper()
		
		for m in regex_compiled.finditer( sequence, overlapping=allow_overlaps ):
			# self, motif, contig, positionStart, strand, regexMatch, molecule='dna'
			tmp = SequenceMotif( motif_seq, contig, start, strand, m, molecule )
			site_list.append( tmp )
	
	return ( input_tuples, site_list )
	
########
# main #
########
if __name__ == "__main__":
	pass
	