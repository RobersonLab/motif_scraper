import os

from motif_scraper import rev_comp, make_degenerate_regex, SequenceMotif, fasta_motif_scan

TEST_FILE = os.path.join( os.path.dirname( __file__ ), 'test.fa' )

def test_fasta_motif_scan_plus():
	############################
	# setup expected responses #
	############################
	plus_strand_response = [ SequenceMotif(  'VCGTACGATS', 'test_a', 0, '+', None, molecule='dna', regex_match_start=1, regex_end_val=10, regex_group_sequence='ACGTACGATC' ), SequenceMotif(  'VCGTACGATS', 'test_a', 0, '+', None, molecule='dna', regex_match_start=231, regex_end_val=240, regex_group_sequence='GCGTACGATG' ) ]
	
	#################
	# run test data #
	#################
	result_plus = fasta_motif_scan( TEST_FILE, ( 'VCGTACGATS', 'test_a', 0, 240, '+' ), regex_ready=False, molecule='dna' )[1]
	
	################
	# assert tests #
	################
	assert( result_plus == plus_strand_response )
	
def test_fasta_motif_scan_minus():
	############################
	# setup expected responses #
	############################
	neg_strand_response = [ SequenceMotif(  'SATCGTACGB', 'test_a', 0, '-', None, molecule='dna', regex_match_start=1, regex_end_val=10, regex_group_sequence='ACGTACGATC' ), SequenceMotif(  'SATCGTACGB', 'test_a', 0, '-', None, molecule='dna', regex_match_start=231, regex_end_val=240, regex_group_sequence='GCGTACGATG' ) ]
	
	#################
	# run test data #
	#################
	result_minus = fasta_motif_scan( TEST_FILE, ( 'SATCGTACGB', 'test_a', 0, 240, '-' ), regex_ready=False, molecule='dna' )[1]
	
	################
	# assert tests #
	################
	assert( result_minus == neg_strand_response )
