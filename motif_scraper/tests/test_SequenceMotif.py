from motif_scraper import rev_comp, make_degenerate_regex, SequenceMotif

def test_SequenceMotif_plus():
	motif_a = SequenceMotif(  'VCGTACGATS', 'test_a', 0, '+', None, molecule='dna', regex_match_start=1, regex_end_val=10, regex_group_sequence='ACGTACGATC' )
	
	assert( str( motif_a ) == "test_a,1,10,+,ACGTACGATC,VCGTACGATS" )

def test_SequenceMotif_minus():
	motif_b = SequenceMotif(  'SATCGTACGB', 'test_a', 0, '-', None, molecule='dna', regex_match_start=1, regex_end_val=10, regex_group_sequence='GATCGTACGT' )
	
	assert( str( motif_b ) == 'test_a,10,1,-,ACGTACGATC,SATCGTACGB' )
