from motif_scraper import rev_comp

def test_rev_comp_dna():
	assert( rev_comp( 'ABCDGHKMNRSTVWY', molecule="dna" ) == 'RWBASYNKMDCHGVT' )

def test_rev_comp_rna():
	assert( rev_comp( 'ABCDGHKMNRSUVWY', molecule="rna" ) == 'RWBASYNKMDCHGVU' )
