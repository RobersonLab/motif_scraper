from motif_scraper import make_degenerate_regex

def test_make_degenerate_regex_dna_all_degen():
	assert( make_degenerate_regex( 'ABCDGHKMNRSTVWY', molecule='dna' ) == 'A[CGT]C[AGT]G[ACT][GT][AC][ACGT][AG][GC]T[ACG][AT][CT]' )

def test_make_degenerate_regex_rna_all_degen():
	assert( make_degenerate_regex( 'ABCDGHKMNRSUVWY', molecule='rna' ) == 'A[CGU]C[AGU]G[ACU][GU][AC][ACGU][AG][GC]U[ACG][AU][CU]' )
	
def test_make_degenerate_regex_dna_repeat_bases():
	assert( make_degenerate_regex( 'ACGTTTTBBNNNNN', molecule='dna' ) == 'ACGT{4}[CGT]{2}[ACGT]{5}' )
